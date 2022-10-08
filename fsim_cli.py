# -*- coding: utf-8 -*-
# pylint: disable=line-too-long
"""
Created on Sat Aug 14 16:50:36 2021

@author: Chris
"""
import json
import argparse
import os
import math
import xml.etree.ElementTree as ET
import subprocess

STD_CFG = os.path.splitext(__file__)[0]+".json"

M = 0.0289644 # molar mass of Earth's air kg/mol
g = 9.80665 # gravitational acceleration: m/s**2
R = 8.3144598 # universal gas constant:  J/(mol·K)

METERS_TO_FEET = 3.28084 # ft/m
FEET_TO_NM = 1.645788e-4 # NM/ft
KNOTS_TO_FPM = 1/FEET_TO_NM/60

def __parse_cg_types(input_str:str):
    try:
        ret = float(input_str)
    except ValueError:
        ret = input_str
    return ret

def _get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--compute-cg",
                        metavar=("ac_icao", "plt_frnt_pssnger", "fuel", "aft_pssngrs", "baggage"),
                        nargs = 5,
                        type = __parse_cg_types,
                        help = "Compute CG, Moment and total mass."
                        )
    parser.add_argument("--slope-gnd",
                        metavar=("DME","AGL","ANGLE"),
                        nargs=3,
                        type = __parse_cg_types,
                        help = "Compute missing parameter for ground based triangle."
                        )
    parser.add_argument("--slope-spd",
                        metavar=("GS","FPM","ANGLE"),
                        nargs=3,
                        type = __parse_cg_types,
                        help = "Compute missing parameter for speed based triangle."
                        )
    parser.add_argument("--slant-corr",
                        metavar=("DME","AGL"),
                        nargs=2,
                        type = float,
                        help = "Calculate corrected distance from DME."
                        )
    parser.add_argument("--dtpp-pdf",
                        metavar=("AIRPORT_ICAO","FILTERS"),
                        nargs="+",
                        help = "Find all d-TPP documents related to airport in a local d-TTP installation."
                        )
    parser.add_argument("--dcs-pdf",
                        metavar=("AIRPORT_ICAO"),
                        nargs="+",
                        help="Find chart supplement given an airport ICAO code."
                        )
    return parser.parse_args()

def _triag(hor:float=None, vert:float=None, angle:float=None):
    """
    Compute missing parameter from a right angled triangle.

    Parameters
    ----------
    hor : float, optional
        Horizontal component, by default None.
    vert : float, optional
        Vertical component, by default None.
    angle : float, optional
        angle in degrees, by default None.

    Returns
    -------
    result : float
        The missing parameter.

    Raises
    ------
    ValueError
        Raised if the number of args is not exactly two.
    """
    if isinstance(hor, (int,float)) and isinstance(vert, (int,float)) and angle is None:
        result = math.atan(vert/hor)*180/math.pi
    elif isinstance(hor, (int,float)) and vert is None and isinstance(angle, (int,float)):
        angle *= math.pi/180 # to radian
        result = math.tan(angle)*hor
    elif hor is None and isinstance(vert, (int,float)) and isinstance(angle, (int,float)):
        angle *= math.pi/180 # to radian
        result = vert*(math.tan(angle))**-1
    else:
        raise ValueError("Function accepts exactly two arguments.")

    return result

def dttp_pdf(args, cfg):
    if len(args.dtpp_pdf) > 1:
        result = get_pdfs_from_dtpp_xml(cfg, args.dtpp_pdf[0], [r.upper() for r in args.dtpp_pdf[1:]])
    else:
        result = get_pdfs_from_dtpp_xml(cfg, args.dtpp_pdf[0])
    for res in result:
        print("[{}] {}: {}, Changed:{}".format(*res))
    if len(args.dtpp_pdf) > 2:
        for i, res in enumerate(result):
                # each time this cmd is run, open a new window
            if (i == 0) and "n" in args.dtpp_pdf[2].lower():
                cmd = "{:s} -new-window {:s}".format(cfg["SUMATRA"], os.path.join(cfg["DTTP_PDF"], res[2]))
            else:
                cmd = "{:s} {:s}".format(cfg["SUMATRA"], os.path.join(cfg["DTTP_PDF"], res[2]))
            subprocess.Popen(cmd, shell=True)


def interp_lin(c1:tuple, c2:tuple, x:float):
    # pylint: disable=invalid-name; shorter names are clearer here
    """
    Performs linear interpolation of a function y(x) between
    to points c1, c2 for a desired point x.

    Parameters
    ----------
    c1 : tuple
        First point (x1,y1) of the linear function.
    c2 : tuple
        Second point (x2,y2) of the linear function.
    x : float
        Independent variable, must be in [x1, x2].

    Returns
    -------
    y : float
        Interpolated value.

    """
    return (c2[1]-c1[1])/(c2[0]-c1[0])*(x-c1[0]) + c1[1]

def conv_vor_hdg(mag_crs:float, mag_dsplcmnt:float, vor_dsplcmnt:float):
    """
    Converts a magnetic course to/from the VOR to a VOR course

    Parameters
    ----------
    mag_crs : float
        Magnetic course to VOR.
    mag_dsplcmnt : float
        Current Magnetic displacment. East -> negative sign.
    vor_dsplcmnt : float
        VOR Magnetic displacment. East -> negative sign.

    Returns
    -------
    mag_vor : float
        Course to/from the VOR.

    """
    return mag_crs + (mag_dsplcmnt - vor_dsplcmnt)

def baro_height(P:float, T_ref:float, P_ref:float, h_ref:float, lapse_rate:float=-6.5e-3):
    """
    Use barometric formula to compute height above MSL.

    Parameters
    ----------
    P : float
        Pressure in mbar to calculate the height at.
    T_ref : float
        Reference Temperature (i.e. METAR) in °C.
    P_ref : float
        Reference Pressure (i.e. METAR) in mbar.
    h_ref : float
        Reference height (i.e. field elevation) in meters.
    lapse_rate : float, optional
        Temperature lapse rate in K/m. The default is -6.5e-3 (ISA).

    Returns
    -------
    h : float
        height where target pressure P has been reached.

    """
    T_ref += 273.15 # to Kelvin
    if math.isclose(lapse_rate, 0, rel_tol=1e-6, abs_tol=1e-8):
        dh = -R*T_ref/(g*M)*math.log(P/P_ref)
    else:
        dh = 1/lapse_rate*( T_ref * (P/P_ref)**(-(R*lapse_rate)/(g*M)) - T_ref)

    return dh + h_ref

def baro_prs(h:float, T_ref:float, P_ref:float, h_ref:float, lapse_rate:float=-6.5e-3):
    """
    Computes barometric pressure at a desired height above ground level.

    Parameters
    ----------
    h : float
        Height to compute the barometric pressure at in meters.
    T_ref : float
        Temperature at reference ground station in °C.
    P_ref : float
        Reference pressure at ground station in millibars/hPa.
    h_ref : float
        Elevation of ground station in meters.
    lapse_rate : float, optional
        Temperature lapse rate in K/m. The default is -6.5e-3 (ISA).

    Returns
    -------
    P : float
        Barometric pressure at desired height in mbar/hPa.

    """

    T_ref += 237.15 # to kelvin
    P_ref *= 100 # to Pascal

    if math.isclose(lapse_rate, 0, rel_tol=1e-6, abs_tol=1e-8):
        P = -g*M*(h-h_ref)/(R*T_ref)
        P = P_ref*math.exp(P)
    else:
        P = ((T_ref+lapse_rate*(h-h_ref))/T_ref)**(-(g*M)/(R*lapse_rate))
        P *= P_ref

    P /= 100 # to millibar / hPa

    return P

def compute_cg(cfg, ac_icao, **kwargs):
    """
    Computes an aircrafts center of gravity (in inches aft of datum),
    total weight (in lbs) and total moment (in lbs-in).

    Parameters
    ----------
    cfg : dict
        Program config
    ac_icao : str
        ICAO identifier of the aircraft. Must be in WEIGHT_CFG.
    **kwargs : float
        Weights of passangers/baggage/fuel. See keys in WEIGHT_CFG.

    Raises
    ------
    KeyError
        Raised if ac_icao is not present in WEIGHT_CFG.

    Returns
    -------
    CG : float
        Center of gravity (in inches aft of datum)
    ac_total_moment : float
        Total moment (in lbs-in)
    ac_total_weight : float
        Total weight (in lbs).
    """
    # select weight dictionary
    if ac_icao.lower() not in cfg["CG_CFGS"].keys():
        raise KeyError("Aircraft {} not implemented.".format(ac_icao))
    mom_dct = cfg["CG_CFGS"][ac_icao]

    # compute moments and total weight:
    ac_total_weight = mom_dct["empty_weight"]
    ac_total_moment = mom_dct["empty_moment"]
    for key, value in kwargs.items():
        if key in mom_dct:
            ac_total_weight += value
            ac_total_moment += interp_lin((0,0), mom_dct[key], value)*1000
        else:
            print("Warning: Keyword {} not found in moment dictionary. Skipping ..".format(key))

    # compute center of gravity
    center_grav = ac_total_moment/ac_total_weight # inches

    return center_grav, ac_total_moment, ac_total_weight

def get_pdfs_from_dcs_xml(cfg:dict, airport_icao:str):
    tree = ET.parse(cfg["DCS_XML"])
    root = tree.getroot()
    if len(airport_icao) == 3:
        cs_filename = root.findall("./location/airport/[aptid='{:s}']/pages/pdf".format(airport_icao.upper()))
    elif len(airport_icao) == 4:
        # omit first character
        cs_filename = root.findall("./location/airport/[aptid='{:s}']/pages/pdf".format(airport_icao[1:].upper()))
    else:
        raise ValueError("Airport ICAO code must have 3 or 4 letters. Code passed was {:s}".format(airport_icao))
    return cs_filename[0].text

def get_pdfs_from_dtpp_xml(cfg:dict, airport_icao:str, code_filters:list=None):
    """
    Parse the xml file provided with the FAA d-TTP publication
    and return all procedure pdfs for a given airport ICAO code.

    Parameters
    ----------
    cfg : dict
        Program config
    airport_icao : str
        Airport ICAO identifier i.e. KAST, KLAX..

    Returns
    -------
    result : list
        List of lists. Each sublist has three result elements:
        Produre type, procedure name, pdf file name
    """
    tree = ET.parse(cfg["DTTP_XML"])
    root = tree.getroot()
    airport_records = root.findall("./state_code/city_name/airport_name/[@icao_ident='{:s}']".format(airport_icao.upper()))
    record_list = []
    for record in airport_records[0]:
        chart_code = record.find("chart_code").text
        if code_filters is None:
            record_list.append([chart_code,
                                record.find("chart_name").text,
                                record.find("pdf_name").text,
                                record.find("cn_flg").text])
        elif chart_code.upper() in code_filters:
            record_list.append([chart_code,
                                record.find("chart_name").text,
                                record.find("pdf_name").text,
                                record.find("cn_flg").text])
        else:
            # doesn't match any filter
            pass
    return record_list

def slant_correction(dme:float, height:float):
    """
    Calculate corrected distance from DME.

    Parameters
    ----------
    dme : float
        DME readout in NM.
    height : float
        Height agl in feet.

    Returns
    -------
    result : float
        Horizontal distance from DME.
    """
    dme /= FEET_TO_NM # to feet
    angle = math.asin(height/dme)*180/math.pi

    return _triag(None, height, angle)*FEET_TO_NM

def slope_gnd(**kwargs):
    """
    Compute Missing parameter from DME, Height AGL and glideslope angle

    Returns
    -------
    result : dict
        Computed results in a dict.
    """
    missing_key = None
    for k,v in kwargs.items():
        if not isinstance(v,(float,int)):
            missing_key = k
            kwargs[k] = None
        elif "dme" in k:
            kwargs[k] /= FEET_TO_NM
    kwargs[missing_key] = _triag(*kwargs.values())
    kwargs["dme"] *= FEET_TO_NM
    return kwargs

def slope_spd(**kwargs):
    """
    Compute Missing parameter from DME, Height AGL and glideslope angle

    Returns
    -------
    result : dict
        Computed results in a dict.
    """
    missing_key = None
    for k,v in kwargs.items():
        if not isinstance(v,(float,int)):
            missing_key = k
            kwargs[k] = None
        elif "gs" in k:
            kwargs[k] *= KNOTS_TO_FPM
    kwargs[missing_key] = _triag(*kwargs.values())
    kwargs["gs"] /= KNOTS_TO_FPM
    return kwargs

def main():
    # pylint: disable=missing-function-docstring
    args = _get_args()
    with open(STD_CFG, "r") as f:
        cfg = json.load(f)

    if args.compute_cg is not None:
        kws = ["ac_icao", "plt_frnt_pssnger", "fuel", "aft_pssngrs", "baggage"]
        kws =  dict(zip(kws, args.compute_cg))
        center_grav, ac_total_moment, ac_total_weight = compute_cg(cfg, **kws)
        print(center_grav, ac_total_moment, ac_total_weight)
    if args.slope_gnd is not None:
        result = slope_gnd(**dict(zip(["dme", "height", "angle"], args.slope_gnd)))
        print("Distance: {:.1f}NM, Height AGL: {:.0f}ft, Angle: {:.1f}deg".format(*result.values()))
    if args.slope_spd is not None:
        result = slope_spd(**dict(zip(["gs", "fpm", "angle"], args.slope_spd)))
        print("Groundspeed: {:.0f}kt, Vert. Speed: {:.0f}ft, Angle: {:.1f}deg".format(*result.values()))
    if args.slant_corr is not None:
        result = slant_correction(*args.slant_corr)
        print("Distance: {:.1f}NM, Height AGL: {:.0f}ft".format(result, args.slant_corr[1]))
    if args.dtpp_pdf is not None:
        dttp_pdf(args, cfg)
    if args.dcs_pdf is not None:
        dcs_pdf = get_pdfs_from_dcs_xml(cfg, args.dcs_pdf[0])
        if len(args.dcs_pdf) > 1:
            if "n" in args.dcs_pdf[1].lower():
                cmd = "{:s} -new-window {:s}".format(cfg["SUMATRA"], os.path.join(cfg["DCS_PDF"], dcs_pdf))
            else:
                cmd = "{:s} {:s}".format(cfg["SUMATRA"], os.path.join(cfg["DCS_PDF"], dcs_pdf))
            subprocess.Popen(cmd, shell=True)
        else:
            print(dcs_pdf)

if __name__ == "__main__":
    main()
