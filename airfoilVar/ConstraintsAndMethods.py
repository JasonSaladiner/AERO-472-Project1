# Constants & Objectives
import numpy as np
import scipy.optimize
import scipy.optimize as sco
from dataclasses import dataclass

cruise_dist_nm = 300  # Nautical miles
cruise_speed_knots = 150  # Knots

loiter_time_hours = 2

cruise_altitude_ft = 12000
cruise_density_slugft3 = 1.6480 * 10 ** -3  # slugs / ft^3

weight_pp_lbf = 1000
minimum_static_margin = 0.15

wing_washout_deg = 2

sea_level_min_speed_knots = 75  # At sea level, weight W0, 60deg F
sea_level_density_slugft3 = 2.3769 * 10 ** -3  # slugs / ft^3

# TAIL FOIL IS NACA0012
extra_fuel_pct = 0.06

req_roc_fpm = 400  # feet per minute
climb_speed_knots = sea_level_min_speed_knots * 1.2
# Must climb at W0.

# Objectives
# 1. Size main wing & tail
#   Main Wing:
#       No sweep
#       2deg washout
#   Tail:
#       Linear or no taper
#       Fixed AR = AR_wing - 2
#       Moving stabilator
# 2. Choose airfoil for main wing
#
# 3. Size main wing plain flap / no flap needed
#

c_f = 0.0045  # Coefficient for viscous drag
C_lbf_per_hrhp = 0.45  # Specific fuel consumption
eta_propeller = 0.85  # Propeller efficiency

fuel_consumption_initial = 0.045  # 4.5% of W0 lbs of fuel required to taxi, takeoff, climbout, etc.
fuel_consumption_landing = 0.01  # 1% of W0 lbs of fuel required to land

eta_tail = 0.9  # Qinf efficiency
kappa = 0.5  # Alpha gradient thingy


def vehicle_len_ft(initial_weight_lbf):
    return 4.47 * (initial_weight_lbf ** 0.23)


def wetted_fuselage_and_vtail_area_ftsquared(vehicle_len_ft):
    return 3.1 * (vehicle_len_ft ** 1.3)


def reymer_deltaclmax(s_ratio, hinge_cos, deltacl_device=0.9):
    """

    :param s_ratio: s_flapped / s_reference
    :param hinge_cos: cos(delta_hinge); hinge angle.
    :param deltacl_device: 0.9 for plain flap
    :return: delta cl max
    """
    return 0.9 * deltacl_device * s_ratio * hinge_cos


def weight_fuselage_lbf(initial_weight_lbf, length_ft):
    return 68.3 * (initial_weight_lbf ** 0.14) * (length_ft ** 0.38)


def weight_wing_lbf(initial_weight_lbf, wing_are_ftsquared, AR_wing):
    return 0.054 * (initial_weight_lbf ** 0.4) * (wing_are_ftsquared ** 0.36) * (AR_wing ** 1.7)


def weight_tail_lbf(initial_weight_lbf, aspect_ratio_tail):
    return 0.028 * (initial_weight_lbf ** 0.9) * (aspect_ratio_tail ** 0.12)


def weight_engine_lbf(max_power_hp):
    return 2.0 * max_power_hp


def weight_flaps(span_ft):
    return 0 if abs(span_ft) < 10 ** -10 else 10 + 3 * span_ft


def x_cg_from_nose_ft(L_fuse_ft, x_wing_quarterchord_ft, w_fuselage_lbf, w_engine_lbf, w_tail_lbf,
                      w_fuel_lbf, w_wing_lbf, w_flaps_lbf):
    w1 = w_fuselage_lbf + w_engine_lbf + w_tail_lbf
    w2 = w_fuel_lbf + w_wing_lbf + w_flaps_lbf
    return ((0.3 * L_fuse_ft * w1) + (x_wing_quarterchord_ft * w2)) / (w1 + w2)


def power_required_hp(weight_lbf, velocity_mph, propeller_efficiency, roc_fps, reference_area_ftsquared, aspect_ratio,
                      oswald_e, dynamic_pressure_psf, CD0):
    v = velocity_mph * 1.46667  # ft/s
    roc = roc_fps  # ft/s
    wovers = weight_lbf / reference_area_ftsquared  # lbf / ft^2
    term1 = (CD0 * dynamic_pressure_psf) / wovers
    term2 = wovers / (np.pi * oswald_e * aspect_ratio * dynamic_pressure_psf)
    term3 = roc / v
    power = (weight_lbf * v) / propeller_efficiency * (term1 + term2 + term3)  # lbf * ft / s
    return power / 550


def range_equation_ft(eta_prop, C_lbf_per_hrhp, LoverD, initial_weight_lbf, final_weight_lbf):
    c = C_lbf_per_hrhp * (1 / 550) * (1 / 3600)  # 1 / ft
    return (eta_prop / c) * LoverD * np.log(initial_weight_lbf / final_weight_lbf)


def fuel_cost_range_lbf(range_miles, eta_prop, C_lbf_per_hrhp, LoverD, initial_weight_lbf):
    c = C_lbf_per_hrhp * (1 / 550) * (1 / 3600)  # 1 / ft
    r = range_miles * 5280  # ft
    exponent = (r * c) / (eta_prop * LoverD)
    w_f = initial_weight_lbf / np.exp(exponent)
    return initial_weight_lbf - w_f


def fuel_cost_endurance_lbf(endurance_hours, eta_prop, loiter_speed_mph, C_lbf_per_hrhp, LoverD, initial_weight_lbf):
    c = C_lbf_per_hrhp * (1 / 550) * (1 / 3600)  # 1 / ft
    v = loiter_speed_mph * 5280  # ft / hr
    E = endurance_hours
    exponent = (E * c * v) / (eta_prop * LoverD)
    w_f = initial_weight_lbf / np.exp(exponent)
    return initial_weight_lbf - w_f


def get_loiter_speed_mph(endurance_hrs, eta_prop, C_lbf_per_hrhp, initial_weight_lbf, min_airspeed_mph, cd,
                     density_sf3, S_ref_ft2):

    c = C_lbf_per_hrhp * (1 / 550) * (1 / 3600)  # 1 / ft
    E = endurance_hrs

    # Note C1, C2, C3 > 0
    def cost_fxn(v):
        v_ftperhr = v * 5280
        cl_req = required_clArea(density_sf3, v*0.868976, initial_weight_lbf) / S_ref_ft2
        LoverD = cl_req / cd(cl_req)
        exponent = (E * c * v_ftperhr) / (eta_prop * LoverD)
        return initial_weight_lbf / np.exp(exponent)

    optimal_v = sco.minimize_scalar(cost_fxn, bounds=[min_airspeed_mph, min_airspeed_mph * 5], method="bounded")
    return optimal_v.x, optimal_v.fun


@dataclass
class FuelBurn:
    total_fuel_weight: float = 0
    taxi_takeoff_fuel: float = 0
    cruise_fuel: float = 0
    loiter_fuel: float = 0
    landing_fuel: float = 0

    w0: float = 0  # Initial
    w1: float = 0  # Post taxi-takeoff
    w2: float = 0  # Post cruise
    w3: float = 0  # Post endurance
    w4: float = 0  # Post landing

    CL_cruise: float = 0
    LoverD_cruise: float = 0
    loiter_speed: float = 0
    CL_loiter: float = 0
    LoverD_loiter: float = 0


def fuel_weight(dry_weight_lbf, cd, sref_ft2, CLMax,
                req_range_miles=cruise_dist_nm * 1.1508,
                eta_prop=eta_propeller,
                endurance_hours=loiter_time_hours,
                min_speed_mph=sea_level_min_speed_knots * 1.15078,
                density_loiter_sf3=cruise_density_slugft3,
                C=C_lbf_per_hrhp,
                ext_res=False):
    burndown = FuelBurn()
    burndown.total_fuel_weight = 0
    while True:
        burndown.w0 = dry_weight_lbf + burndown.total_fuel_weight
        # 0 -> 1: taxi, takeoff, ...
        burndown.taxi_takeoff_fuel = fuel_consumption_initial * burndown.w0
        burndown.w1 = burndown.w0 - burndown.taxi_takeoff_fuel
        # 1 -> 2: range
        burndown.CL_cruise = required_clArea(cruise_density_slugft3, cruise_speed_knots, burndown.w1) / sref_ft2
        burndown.LoverD_cruise = burndown.CL_cruise / cd(burndown.CL_cruise)
        burndown.cruise_fuel = fuel_cost_range_lbf(req_range_miles, eta_prop, C, burndown.LoverD_cruise, burndown.w1)
        burndown.w2 = burndown.w1 - burndown.cruise_fuel
        # 2 -> 3: endurance
        # ls_mph, fuel_23_lbf = get_loiter_speed_mph(endurance_hours, eta_prop, C, w2, min_speed_mph, cd,
        #                                            density_loiter_sf3, sref_ft2)

        ls_mph = required_airspeed_mph(density_loiter_sf3, CLMax * sref_ft2, burndown.w2)

        def min_fuel_fxn(v_mph):
            cl_k = required_clArea(density_loiter_sf3, v_mph*0.868976, burndown.w2) / sref_ft2
            # print("V: {} mph -> CL: {}".format(v_mph, cl_k))
            ld_k = cl_k / cd(cl_k)
            return fuel_cost_endurance_lbf(endurance_hours, eta_prop, v_mph, C, ld_k, burndown.w2)

        burndown.loiter_speed = scipy.optimize.minimize_scalar(min_fuel_fxn, bounds=[ls_mph, ls_mph * 5],
                                                               method="bounded").x

        burndown.CL_loiter = required_clArea(density_loiter_sf3, burndown.loiter_speed*0.868976, burndown.w2) / sref_ft2
        burndown.LoverD_loiter = burndown.CL_loiter / cd(burndown.CL_loiter)
        burndown.loiter_fuel = fuel_cost_endurance_lbf(endurance_hours, eta_prop, burndown.loiter_speed,
                                                       C, burndown.LoverD_loiter, burndown.w2)

        burndown.w3 = burndown.w2 - burndown.loiter_fuel
        # 3 -> 4: landing
        burndown.landing_fuel = fuel_consumption_landing * burndown.w0
        burndown.w4 = burndown.w3 - burndown.landing_fuel
        # 4: 6% fuel reserve
        fuel_extra = burndown.total_fuel_weight * extra_fuel_pct
        fuel_actual = burndown.taxi_takeoff_fuel + burndown.cruise_fuel + burndown.loiter_fuel +\
                      burndown.landing_fuel + fuel_extra
        offset = burndown.total_fuel_weight - fuel_actual
        if abs(offset) < 10 ** -5:
            if ext_res:
                return burndown
            return burndown.total_fuel_weight, burndown.loiter_speed
        burndown.total_fuel_weight -= offset


def required_clArea(rho_slugft3, v_min_knots, lift_lbf):
    # v = v_min_knots * 1.68781
    # q = 0.5 * rho_slugft3 * (v ** 2)
    # return lift_lbf / q
    rho = 515.378819 * rho_slugft3
    v = v_min_knots * 0.514444
    q = 0.5 * rho * v * v
    F = lift_lbf * 4.44822
    return (F / q) * 10.7639


def required_airspeed_mph(rho_slugft3, CLArea, lift_lbf):
    # v = np.sqrt(2 * lift_lbf / (rho_slugft3 * CLArea))
    # return v
    F = lift_lbf * 4.44822
    rho = 515.378819 * rho_slugft3
    cla = CLArea / 10.7639
    return np.sqrt((2 * F) / (rho * cla)) * 2.23694



def required_power(initial_weight_lbf, wing_area_ft2, AR, oswald_e, CD0):
    v = climb_speed_knots * 1.68781  # ft/s
    v_mph = v * 0.681818
    q = 0.5 * sea_level_density_slugft3 * (v ** 2)
    cr = req_roc_fpm / 60  # ft / s
    return power_required_hp(initial_weight_lbf, v_mph, eta_propeller, cr, wing_area_ft2, AR, oswald_e, q, CD0)


def aerodynamic_center_approx(eta_tail, tail_len, a_tail, kappa, tail_area, fuse_moment_coeff,
                              a_wing, wing_area):
    """

    :param eta_tail: qinf efficiency
    :param tail_len: dimensional, ft, wrt main wing AC
    :param a_tail: a = CL_alpha  S derivative, 1/radians
    :param kappa: down wash gradient
    :param tail_area: ft ^ 2
    :param fuse_moment_coeff: swirly A in ppt, A = 0.6 d^2 *lf ^2 +- 50%; mf = 2qA(alpha)
    :param a_wing: 1/radians
    :param wing_area: ft ^ 2
    :return: x ac (ft), w.r.t main wing AC
    """

    numerator = (eta_tail * a_tail * (1 - kappa) * tail_area * tail_len) - (2 * fuse_moment_coeff)
    dnom = (a_wing * wing_area) + (eta_tail * a_tail * (1 - kappa) * tail_area)
    return numerator / dnom


def compute_main_wing_pos_from_nose(static_margin, l_fuse_ft, w_fuselage_engine_tail_lbf, w_wing_flap_lbf, w_fuel_lbf,
                                    chord_ft, eta_tail, a_tail, kappa, fuse_mom_coeff, a_wing,
                                    wing_area, tail_volume_coeff, full=False):
    x_qc = 0
    while True:
        tail_arm_len_ft = l_fuse_ft - x_qc
        tail_area_ft2 = wing_area * chord_ft * tail_volume_coeff / tail_arm_len_ft
        x_cg_k_from_nose = x_cg_from_nose_ft(l_fuse_ft, x_qc, w_fuselage_engine_tail_lbf, 0, 0, w_fuel_lbf,
                                             w_wing_flap_lbf, 0)
        ac_k_from_nose = aerodynamic_center_approx(eta_tail, tail_arm_len_ft, a_tail, kappa, tail_area_ft2,
                                                   fuse_mom_coeff, a_wing, wing_area) + x_qc
        sm_k = (ac_k_from_nose - x_cg_k_from_nose) / chord_ft
        delta_sm = sm_k - static_margin
        if abs(delta_sm) < 10 ** -10:
            if full:
                return x_qc, tail_area_ft2, x_cg_k_from_nose, ac_k_from_nose
            return x_qc, tail_area_ft2
        x_qc -= (delta_sm * chord_ft)


def cd0(c_f, cdi_0, s_wet, s_ref):
    cd_v = s_wet / s_ref * c_f
    return cd_v + cdi_0



if __name__ == "__main__":


    dry_weight = 3000
    loverd = 10
    aspect_ratio = 8

    w_fuel_lbf = fuel_weight(dry_weight, loverd, lambda x: 0.03, 30, 0.8, aspect_ratio)[0]
    print("Fuel weight for {} lbf plane w/ max L/D of {}: {}".format(dry_weight, loverd, w_fuel_lbf))

    total_weight_lbf = dry_weight + w_fuel_lbf
    CLA = required_clArea(sea_level_density_slugft3, sea_level_min_speed_knots, total_weight_lbf)
    print("Required CL*A: {} ft^2".format(CLA))

    oswald_e = 0.8
    area_ft2 = CLA
    CD0 = 0.03
    target_tvc = 0.5
    span_ft = np.sqrt(aspect_ratio * area_ft2)
    chord_ft = area_ft2 / span_ft

    p_req = required_power(total_weight_lbf, CLA, aspect_ratio, oswald_e, CD0)
    print("Required power: {} hp".format(p_req))

    w_wing = weight_wing_lbf(total_weight_lbf, area_ft2, aspect_ratio)
    w_flaps = 0
    w_engine = weight_engine_lbf(p_req)
    w_tail = weight_tail_lbf(total_weight_lbf, aspect_ratio - 2)
    fuse_len_ft = vehicle_len_ft(total_weight_lbf)
    w_fuse = weight_fuselage_lbf(total_weight_lbf, fuse_len_ft)

    w_dry = weight_pp_lbf + w_wing + w_flaps + w_engine + w_tail + w_fuse

    print("Dry Weight buildup: PP + Engine + Wing + Tail + Fuselage = {} + {} + {} + {} + {} = {} lbf".format(
        weight_pp_lbf, w_wing, w_engine, w_tail, w_fuse, w_dry
    ))

    x_mw_ac, tail_area_ft2 = compute_main_wing_pos_from_nose(minimum_static_margin, fuse_len_ft,
                                                             w_engine + w_fuse + w_tail, w_wing + w_flaps,
                                                             w_fuel_lbf, chord_ft,
                                                             eta_tail,
                                                             np.pi * 2,  # HS cl alpha
                                                             kappa,
                                                             0,
                                                             np.pi * 2,  # MW cl alpha
                                                             area_ft2, target_tvc)
    print("Fuse len: {} ft".format(fuse_len_ft))
    print("Main wing quarter cord position: {} ft, HS area {} ft^2".format(x_mw_ac, tail_area_ft2))

    results = []
    settings = []

    loiter_speed_mph = -1
    cl_max_foil = 1
    cdi_0 = 0
    for flap_pct in [0, 0.1, 0.2, 0.3, 0.4]:
        for aspect_ratio in [5, 6, 7, 8, 9, 10, 11, 12]:
            print("ASPECT RATIO: {}".format(aspect_ratio))
            oswald_e = 0.8  # TODO: find values

            # Dry weight loop
            loverd_cruise = None  # Modified in loop; exists to prevent printing error
            dry_weight = 3000
            dry_weight_km1 = 0.0
            area_ft2 = required_clArea(sea_level_density_slugft3, sea_level_min_speed_knots, dry_weight) / cl_max_foil
            CD0 = 0.03
            while abs(dry_weight - dry_weight_km1) > 10 ** -5:
                # FUEL WEIGHT CALCS
                # TODO: tail aerodynamic effects in LoverD
                CL_cruise = required_clArea(cruise_density_slugft3, cruise_speed_knots, total_weight_lbf) / area_ft2
                loverd_cruise = CL_cruise / (CD0 + ((CL_cruise ** 2) / (np.pi * oswald_e * aspect_ratio)))

                def CD_cl(cl):
                    return CD0 + ((cl ** 2) / (np.pi * oswald_e * aspect_ratio))

                w_fuel_lbf, loiter_speed_mph = fuel_weight(dry_weight, loverd_cruise, CD_cl, area_ft2)
                total_weight_lbf = dry_weight + w_fuel_lbf

                # MIN SPEED CALCS
                CLA = required_clArea(sea_level_density_slugft3, sea_level_min_speed_knots, total_weight_lbf)
                delta_cl_max = reymer_deltaclmax(flap_pct, 1)

                # ASSUMPTION: elliptical wing -> elliptical cl alpha correction
                a_mw = (np.pi * 2 * aspect_ratio) / (aspect_ratio + 2)

                # ASSUMPTION: hs CLa is 2pi, ar = mw ar - 2 (prescribed)
                hs_ar = aspect_ratio - 2
                a_hs = (np.pi * 2 * hs_ar) / (hs_ar + 2)
                cl_max = cl_max_foil + delta_cl_max
                # TODO: best area assumed to be minimum area
                area_ft2 = CLA / cl_max

                target_tvc = 0.5  # TODO: Arbitrary?
                span_ft = np.sqrt(aspect_ratio * area_ft2)
                chord_ft = area_ft2 / span_ft
                flap_span = span_ft * flap_pct

                # Required Power Calculations
                p_req = required_power(total_weight_lbf, CLA, aspect_ratio, oswald_e, CD0)

                # Weight approximations
                w_wing = weight_wing_lbf(total_weight_lbf, area_ft2, aspect_ratio)
                w_flaps = weight_flaps(flap_span)
                w_engine = weight_engine_lbf(p_req)
                w_tail = weight_tail_lbf(total_weight_lbf, aspect_ratio - 2)
                fuse_len_ft = vehicle_len_ft(total_weight_lbf)
                w_fuse = weight_fuselage_lbf(total_weight_lbf, fuse_len_ft)
                w_dry = weight_pp_lbf + w_wing + w_engine + w_tail + w_fuse + w_flaps

                # AC and Tail Area calcs (AKA Longitudinal stability requirements)
                x_mw_ac, tail_area_ft2 = compute_main_wing_pos_from_nose(minimum_static_margin, fuse_len_ft,
                                                                         w_engine + w_fuse + w_tail,
                                                                         w_wing + w_flaps,
                                                                         w_fuel_lbf,
                                                                         chord_ft,
                                                                         eta_tail,
                                                                         a_hs,  # HS cl alpha
                                                                         kappa,
                                                                         0,
                                                                         a_mw,  # MW cl alpha
                                                                         area_ft2, target_tvc)

                # correct aerodynamic performance variables from longitudinal calcs (i.e.; L/D, CD0)
                s_wet_ft2 = wetted_fuselage_and_vtail_area_ftsquared(fuse_len_ft) + (2 * area_ft2)\
                            + (2 * tail_area_ft2)

                CD0 = cd0(c_f, cdi_0, s_wet_ft2, area_ft2)

                # Loop variables
                dry_weight_km1 = dry_weight
                dry_weight = w_dry

            print("Required CL*A: {} ft^2, CD0: {}".format(CLA, CD0))
            print("Required power: {} hp".format(p_req))
            print("Fuse len: {} ft".format(fuse_len_ft))
            print("Dry Weight buildup: PP + Engine + Wing + Flaps + Tail + Fuselage = {} + {} + {} + {} + {} + {} = "\
                  "{} lbf".format(
                weight_pp_lbf, w_wing, w_flaps, w_engine, w_tail, w_fuse, w_dry
            ))
            print("Main wing quarter cord position: {} ft, HS area {} ft^2".format(x_mw_ac, tail_area_ft2))
            print("Optimal loiter speed: {} mph".format(loiter_speed_mph))
            print("Fuel weight for {} lbf plane w/ max L/D of {}: {}".format(dry_weight, loverd_cruise, w_fuel_lbf))
            results.append(w_fuel_lbf)
            settings.append([flap_pct, aspect_ratio])

    print("Results:")
    best_key = []
    best_val = 9999999999
    for key, val in zip(settings, results):
        print("{}% flap, {} AR: {} lbf fuel".format(key[0] * 100, key[1], val))
        if val < best_val:
            best_val = val
            best_key = key
    print("BEST: {}% flap, {} AR: {} lbf fuel".format(best_key[0] * 100, best_key[1], best_val))
