import numpy as np

from ConstraintsAndMethods import *

def interpolate(x, x_list, y_list):
    for i in range(1, len(x_list)):
        if x_list[i - 1] < x <= x_list[i]:
            return (y_list[i] - y_list[i - 1]) / (x_list[i] - x_list[i - 1]) * (x - x_list[i - 1]) + y_list[i - 1]
    if x > max(x_list):
        return (y_list[-1] - y_list[-2]) / (x_list[-1] - x_list[-2]) * (x - x_list[-2]) + y_list[-2]
    return None


class WingFoil:

    def __init__(self):
        self.cd0 = 0.01
        self.cl_max = 1

    def cd(self, cl):
        return self.cd0

    def ClMax(self):
        return self.cl_max


class Wing:

    def __init__(self):
        self.AR = 4
        self.e = 0.8
        self.area = 10
        self.airfoil = WingFoil()

    def clalpha(self):
        return (np.pi * 2 * self.AR) / (self.AR + 2)

    def span(self):
        return np.sqrt(self.AR * self.area)

    def mac(self):
        return self.area / self.span()

    def wetted_area(self):
        return self.area * 2

    def ac_le(self):
        return self.mac() * 0.25

    def CLmax(self):
        return self.airfoil.ClMax()

    def CD_noskinfric(self, CL):
        return self.airfoil.cd(CL) + (CL ** 2) / (np.pi * self.e * self.AR)


class Plane:

    def __init__(self):
        # Requirements
        self.static_margin = minimum_static_margin
        self.eta_tail = eta_tail
        self.kappa_tail = kappa
        self.payload = weight_pp_lbf
        self.visc_coeff_fric = c_f

        # Independents
        self.mw = Wing()
        self.hs = Wing()
        self.flap_pct = 0
        self.tail_vol_coeff = 0.5

        # Computed
        self.power_req = 0
        self.fuse_len = 10

        self.wetted_area = 0
        self.wetted_area_nowings = 0
        self.oswald_e = 0.8

        self.fuse_weight = 0
        self.wing_weight = 0
        self.tail_weight = 0
        self.engine_weight = 0
        self.flap_weight = 0

        self.x_ac_ft = 0
        self.x_cg_ft = 0

        self.mw_x_ac = 10
        self.hs_x_ac = 20

    def cd0(self):
        cd_v = self.wetted_area / self.mw.area * self.visc_coeff_fric
        return cd_v + self.mw.airfoil.cd(0) + (self.hs.airfoil.cd(0) * self.eta_tail * (self.hs.area / self.mw.area))

    def dry_weight(self):
        return self.wing_weight + self.flap_weight + self.engine_weight + self.fuse_weight + self.tail_weight \
               + self.payload

    def CLmax(self):
        delta_cl_max = reymer_deltaclmax(self.flap_pct, 1)
        cl_max_mw = self.mw.CLmax() + delta_cl_max
        arm_ratio = (self.mw_x_ac - self.x_cg_ft) / (self.hs_x_ac - self.x_cg_ft)
        return cl_max_mw * (1 - arm_ratio)

    def CD(self, CL):
        arm_ratio = (self.mw_x_ac - self.x_cg_ft) / (self.hs_x_ac - self.x_cg_ft)
        CL_hs = (arm_ratio * CL) * (1 / self.eta_tail) * (self.mw.area / self.hs.area)
        cd_visc = self.wetted_area / self.mw.area * self.visc_coeff_fric
        cd_tail = self.eta_tail * (self.hs.area / self.mw.area) * self.hs.CD_noskinfric(CL_hs)
        return self.mw.CD_noskinfric(CL) + cd_tail + cd_visc

    def build_optimal_airplane(self, fuel_weight):
        total_weight_lbf = self.dry_weight() + fuel_weight
        if total_weight_lbf == 0:
            total_weight_lbf = 3000
        last_cl_max = None
        self.hs.AR = self.mw.AR - 2
        cl_max_k = self.CLmax()
        while last_cl_max is None or abs(cl_max_k - last_cl_max) > 10 ** -3:
            CLA = required_clArea(sea_level_density_slugft3, sea_level_min_speed_knots, total_weight_lbf)

            self.mw.area = CLA / cl_max_k
            mw_span = self.mw.span()
            flap_span = self.flap_pct * mw_span

            self.power_req = required_power(total_weight_lbf, self.mw.area, self.mw.AR, self.oswald_e, self.cd0())
            self.fuse_len = vehicle_len_ft(total_weight_lbf)

            self.wing_weight = weight_wing_lbf(total_weight_lbf, self.mw.area, self.mw.AR)
            self.flap_weight = weight_flaps(flap_span)
            self.engine_weight = weight_engine_lbf(self.power_req)
            self.tail_weight = weight_tail_lbf(total_weight_lbf, self.hs.AR)
            self.fuse_weight = weight_fuselage_lbf(total_weight_lbf, self.fuse_len)

            # Static margin calculations
            self.mw_x_ac, self.hs.area, self.x_cg_ft, \
            self.x_ac_ft = compute_main_wing_pos_from_nose(self.static_margin, self.fuse_len,
                                                           self.engine_weight + self.fuse_weight + self.tail_weight,
                                                           self.wing_weight + self.flap_weight,
                                                           fuel_weight,
                                                           self.mw.mac(),
                                                           self.eta_tail, self.hs.clalpha(), self.kappa_tail,
                                                           0,
                                                           self.mw.clalpha(),
                                                           self.mw.area, self.tail_vol_coeff, full=True)
            self.hs_x_ac = self.fuse_len  # Assumption; see SM calcs
            self.wetted_area = wetted_fuselage_and_vtail_area_ftsquared(self.fuse_len) + self.mw.wetted_area() \
                               + self.hs.wetted_area()
            last_cl_max = cl_max_k
            cl_max_k = self.CLmax()

    def __str__(self):
        return """Weight: {} lbf
    Fuselage: {} lbf
    Wing: {} lbf 
    Flaps: {} lbf ({}%)
    Tail: {} lbf
    Engine: {} lbf
    Payload: {} lbf
CL Max: {} (L/D: {}) 
Power Req: {} hp
X CG: {} ft, X AC: {} ft (from nose) 
    X MW AC: {} ft 
    X HS AC: {} ft
MW Area: {} ft^2 (MAC: {} ft)
    MW AR: {}
    MW Span: {} ft 
HS Area: {} ft^2
    HS AR: {}
    HS Span: {} ft 
        """.format(
            self.dry_weight(),
            self.fuse_weight, self.wing_weight,
            self.flap_weight, self.flap_pct * 100,
            self.tail_weight, self.engine_weight, self.payload,
            self.CLmax(), self.CLmax() / self.CD(self.CLmax()),
            self.power_req,
            self.x_cg_ft, self.x_ac_ft, self.mw_x_ac, self.hs_x_ac,
            self.mw.area, self.mw.mac(), self.mw.AR, self.mw.span(),
            self.hs.area, self.hs.AR, self.hs.span()
        )


naca_2412_cl = [-0.855,
                -0.8389,
                -0.8199,
                -0.8012,
                -0.7832,
                -0.7655,
                -0.7471,
                -0.728,
                -0.708,
                -0.6871,
                -0.6705,
                -0.6492,
                -0.6271,
                -0.5952,
                -0.5616,
                -0.5285,
                -0.4933,
                -0.4591,
                -0.4258,
                -0.3914,
                -0.3562,
                -0.3216,
                -0.2881,
                -0.2559,
                -0.2227,
                -0.1916,
                -0.1628,
                -0.1331,
                -0.106,
                -0.0784,
                -0.0517,
                -0.0252,
                0.0015,
                0.0283,
                0.0551,
                0.0819,
                0.1087,
                0.1353,
                0.1615,
                0.1875,
                0.2123,
                0.2341,
                0.2592,
                0.2995,
                0.3336,
                0.3593,
                0.3851,
                0.4108,
                0.4368,
                0.4629,
                0.4888,
                0.5148,
                0.5404,
                0.5646,
                0.5883,
                0.6121,
                0.637,
                0.6631,
                0.6898,
                0.7165,
                0.7431,
                0.7695,
                0.7956,
                0.8213,
                0.8465,
                0.8697,
                0.8961,
                0.9227,
                0.9491,
                0.9752,
                1.0008,
                1.0235,
                1.0508,
                1.0773,
                1.0995,
                1.1243,
                1.1483,
                1.1718,
                1.1912,
                1.2107,
                1.2331,
                1.2545,
                1.2753,
                1.2952,
                1.3133,
                1.3251,
                1.3301,
                1.3452,
                1.3576,
                1.3658,
                1.3699,
                1.3719,
                ]

naca_2412_cd = [0.0266,
                0.02485,
                0.02379,
                0.02264,
                0.0213,
                0.01988,
                0.01859,
                0.01749,
                0.0166,
                0.01595,
                0.01379,
                0.01309,
                0.01262,
                0.01217,
                0.01178,
                0.0108,
                0.01035,
                0.01008,
                0.00958,
                0.00944,
                0.00942,
                0.00891,
                0.00874,
                0.00859,
                0.00847,
                0.008,
                0.00774,
                0.00752,
                0.00734,
                0.00715,
                0.00686,
                0.00658,
                0.00635,
                0.00615,
                0.00597,
                0.00579,
                0.0056,
                0.00534,
                0.00502,
                0.00475,
                0.00444,
                0.00419,
                0.00417,
                0.0043,
                0.00448,
                0.00464,
                0.00481,
                0.00501,
                0.00521,
                0.00541,
                0.00566,
                0.00593,
                0.00625,
                0.0068,
                0.00746,
                0.00813,
                0.00865,
                0.00897,
                0.00919,
                0.00942,
                0.00965,
                0.00992,
                0.01023,
                0.01059,
                0.01104,
                0.01179,
                0.01199,
                0.01216,
                0.01235,
                0.01258,
                0.01286,
                0.01358,
                0.0136,
                0.01372,
                0.01449,
                0.01484,
                0.01527,
                0.01576,
                0.01681,
                0.0178,
                0.01836,
                0.01902,
                0.01972,
                0.0205,
                0.02144,
                0.02309,
                0.02552,
                0.02657,
                0.02781,
                0.02918,
                0.0307,
                0.03238]

if __name__ == "__main__":

    naca_2412 = WingFoil()
    naca_2412.cd = lambda cl_k: interpolate(cl_k, naca_2412_cl, naca_2412_cd)
    naca_2412.cl_max = 1

    best_plane = None
    best_fuel_weight = None
    best_ldc = None
    best_lde = None
    for fpct in [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        print("FLAP: {}".format(fpct))
        for ar in [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]:
            plane_0 = Plane()
            plane_0.flap_pct = fpct
            plane_0.mw.AR = ar
            plane_0.build_optimal_airplane(0)
            plane_0.mw.airfoil = naca_2412
            plane_0.hs.airfoil = naca_2412
            prev_dry_weight = None
            w_fuel_lbf = 0
            while prev_dry_weight is None or abs(prev_dry_weight - plane_0.dry_weight()) > 10 ** -5:
                prev_dry_weight = plane_0.dry_weight()
                w_fuel_lbf, loiter_speed_mph, ldc, lde = fuel_weight(plane_0.dry_weight(), plane_0.CD,
                                                                     plane_0.mw.area, plane_0.CLmax(), ext_res=True)
                plane_0.build_optimal_airplane(w_fuel_lbf)
            if best_fuel_weight is None or best_fuel_weight > w_fuel_lbf:
                best_fuel_weight = w_fuel_lbf
                best_plane = plane_0
                best_ldc = ldc
                best_lde = lde
    print(best_plane)
    print("FUEL WEIGHT: {}".format(best_fuel_weight))
    print("L/D Cruise: {}, L/D Endurance {}".format(best_ldc, best_lde))
