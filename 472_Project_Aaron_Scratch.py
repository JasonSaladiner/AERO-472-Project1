import numpy as np
import pandas as pd

# slug/ft^3
SL_density = .0024
alt_density = .0016
# ft lbf/sec
hp = 550
# ft
nm = 6076

# totally made up, is actually dependent on We, Wpp, Wf 3168
Wo = 2131.35
Wpp = 1000


def We(Wo):
    val = Wo * 2.36 * Wo ** (-.18)
    return val


def Wfo(Wo, isRefined):
    if isRefined:
        tail_dimensions, wing_dimensions = wingAndTailSizing()[0:5], wingAndTailSizing()[5:]

        P_max, L = power(Wo, minStallArea(Wo), wing_dimensions[2])
        W_empty = refinedEmptyWeight(Wo, minStallArea(Wo), wing_dimensions[2], tail_dimensions[4], P_max, L)
    elif not isRefined:
        W_empty = We(Wo)
    val = Wo - W_empty - Wpp
    return val


def minStallArea(Wo):
    # stall speed check
    V_stall = 75 * 1.6878
    q_stall = .5 * SL_density * V_stall ** 2
    CL_max = 1.2  # no flaps, guessed an attainable value
    S_wing = Wo / (q_stall * CL_max)
    # print(S_wing)
    return S_wing


def initialWeighing(Wo, isRefined):
    if isRefined:
        tail_dimensions, wing_dimensions = wingAndTailSizing()[0:5], wingAndTailSizing()[5:]

        P_max, L = power(Wo, minStallArea(Wo), wing_dimensions[2])
        W_empty = refinedEmptyWeight(Wo, minStallArea(Wo), wing_dimensions[2], tail_dimensions[4], P_max, L)
    elif not isRefined:
        W_empty = We(Wo)
    # climbout from station 0 --> 1 (everything up to cruise phase)
    W1 = .955 * Wo  # after warmup, taxi, etc.
    Wf1 = W1 - Wpp - W_empty  # Weight of fuel at start of cruise
    # print(Wf1)

    # cruise to loiter 1 --> 2
    R = 300 * nm  # nm to ft
    C = 0.45 / hp / 3600  # 1/ft gross units
    n_prop = .85
    L_D_Cruise = 12
    V_cruise = 150*1.6878
    S_wing = minStallArea(Wo)
    q_cruise = .5*alt_density*V_cruise**2
    CL_cruise1 = W1/(q_cruise*S_wing)
    CD_cruise1 = CL_cruise1/L_D_Cruise
    #print(CL_cruise1)
    x = (R * C) / (n_prop * L_D_Cruise)
    W2 = W1 / np.exp(x)  # Weight at end of cruise phase
    CL_cruise2 = W2/(q_cruise*S_wing)
    CD_cruise2 = CL_cruise2/L_D_Cruise
    #print(CL_cruise2)
    Wf2 = W2 - Wpp - W_empty  # Weight of fuel left at end of cruise
    # print(Wf2)

    # loiter 2 --> 3
    E = 7200  # seconds, is endurance
    L_D_loiter = 10
    CL_loiter1 = .8  # guess
    CD_loiter1 = CL_loiter1/L_D_loiter
    q_loiter = W2 / (S_wing * CL_loiter1)
    #print(cruise())
    V_loiter = np.sqrt(2 * q_loiter / alt_density)
    x = E * C * V_loiter / (n_prop * L_D_loiter)
    W3 = W2 / np.exp(x)  # weight at end of loiter
    CL_loiter2 = W3/(q_loiter*S_wing)
    CD_loiter2 = CL_loiter2/L_D_loiter
    #print(CL_loiter2)
    Wf3 = W3 - Wpp - W_empty
    # print(Wf3) # psf

    # descent, landing, tax 3 --> 4
    W4 = .99 * W3
    #print(loiter())
    Wf4 = W4 - Wpp - W_empty
    final_fuel = 100*Wf4 / Wfo(Wo, isRefined)
    return W1, W2, W3, Wf4, final_fuel, CL_cruise1, CL_cruise2, CD_cruise1, CD_cruise2, CL_loiter1, CL_loiter2, CD_loiter1, CD_loiter2


def wingAndTailSizing():
    CHT = .7
    S_wing = minStallArea(Wo)
    AR_wing = 6 # just a good guess
    b_wing = np.sqrt(AR_wing*S_wing)
    c_bar = S_wing/b_wing # assuming a rectangular wing here
    l_tail = .5*b_wing # totally made up relationship, just a guess
    S_tail = CHT*c_bar*S_wing/l_tail
    #print(S_tail)
    St_Sw = S_tail/S_wing
    AR_tail = AR_wing - 2
    b_tail = np.sqrt(AR_tail*S_tail)
    c_tail = S_tail/b_tail
    #print(b_tail)
    return S_tail, l_tail, b_tail, c_tail, AR_tail, b_wing, c_bar, AR_wing


def power(Wo, S_wing, AR_wing):
    e = .85
    RoC = 400 / 60
    V_stall = 75 * 1.6878
    n_prop = .85
    # need to get Cdo
    C_f = .0045
    t_c = 1
    L = 4.47*Wo**.23
    #L = S_wing/wingAndTailSizing()[5]
    S_exposed = .8 * S_wing
    S_wet = 2*(1 + .2*t_c)*S_exposed + 3.1*L**1.3
    Cdo = C_f*S_wet/S_wing
    V_climb = 1.2*V_stall
    q_climb = .5*SL_density*V_climb**2
    P_max = (Wo*V_climb/n_prop)*(Cdo*q_climb/(Wo/S_wing) + (Wo/S_wing)/(np.pi*e*AR_wing*q_climb) + RoC/V_climb)/(hp*32.1741) # gonna want to check this, just kinda guessed the hp term should have a gravity next to it
    #print(P_max)
    return P_max, L


def refinedEmptyWeight(Wo, S_wing, AR_wing, AR_tail, P_max, L):
    L_flaps = 0
    W_fuse = 68.3 * (Wo**.14) * (L**.38)
    #W_fuse = 0
    W_wing = .054 * (Wo**.4) * (S_wing**.36) * (AR_wing**1.7)
    W_tail = .028 * (Wo**.9) * (AR_tail**.12)
    W_eng = 2.0*P_max
    # if flaps are included specify the total length they occupy
    W_flaps = 0*(30 + 3*L_flaps)
    We = W_fuse + W_wing + W_tail + W_eng + W_flaps
    return We


isRefined = False
end_fuel_weight, fuel_percent, cl_cruise_start, cl_cruise_end, cd_cruise_start, cd_cruise_end, cl_loiter_start, cl_loiter_end, cd_loiter_start, cd_loiter_end = initialWeighing(Wo, isRefined)[3:]
tail_dimensions, wing_dimensions = wingAndTailSizing()[0:5], wingAndTailSizing()[5:]
P_max, L = power(Wo, minStallArea(Wo), wing_dimensions[2])
We_revised = refinedEmptyWeight(Wo, minStallArea(Wo), wing_dimensions[2], tail_dimensions[4], P_max, L)
df1 = pd.DataFrame([cl_cruise_start, cl_cruise_end, cd_cruise_start, cd_cruise_end, cl_loiter_start, cl_loiter_end, cd_loiter_start, cd_loiter_end])
df1.index = ['cl_cruise_start', 'cl_cruise_end', 'cd_cruise_start', 'cd_cruise_end', 'cl_loiter_start', 'cl_loiter_end', 'cd_loiter_start', 'cd_loiter_end']
print(f'Pre-revision:\n'
      f'{df1.to_string(header=False)} \n')


# not changing old parameters but updating only the fuel related terms and carrying through
isRefined = True
refined_end_fuel, refined_fuel_percent, cl_cruise_start, cl_cruise_end, cd_cruise_start, cd_cruise_end, cl_loiter_start, cl_loiter_end, cd_loiter_start, cd_loiter_end = initialWeighing(Wo, isRefined)[3:]
df2 = pd.DataFrame([Wo, We(Wo), Wpp, Wfo(Wo, isRefined=False), end_fuel_weight, fuel_percent, minStallArea(Wo), wing_dimensions[0], wing_dimensions[1], wing_dimensions[2], tail_dimensions[0], tail_dimensions[1], tail_dimensions[2], tail_dimensions[3], tail_dimensions[4], We_revised, Wfo(Wo, isRefined), refined_end_fuel, refined_fuel_percent])
df2.index = ['Wo', 'We', 'Wpp', 'Wfo', 'Wf4', 'Percent fuel reserve', 'S_wing', 'b_wing', 'c_bar', 'AR_wing', 'S_tail', 'l_tail', 'b_tail', 'c_tail', 'AR_tail', 'Revised We', 'Refined Wfo', 'Refined Wf4', 'Refined percent fuel reserve']

df3 = pd.DataFrame([cl_cruise_start, cl_cruise_end, cd_cruise_start, cd_cruise_end, cl_loiter_start, cl_loiter_end, cd_loiter_start, cd_loiter_end])
df3.index = ['cl_cruise_start', 'cl_cruise_end', 'cd_cruise_start', 'cd_cruise_end', 'cl_loiter_start', 'cl_loiter_end', 'cd_loiter_start', 'cd_loiter_end']
print(f'Post-revision:\n'
      f'{df3.to_string(header=False)} \n')

# the output suggests the initial sizing is way too large
print(f'Complete table of values:\n'
      f'{df2.to_string(header=False)}')
#print(We_revised)

