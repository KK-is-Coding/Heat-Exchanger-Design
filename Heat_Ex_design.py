import CoolProp.CoolProp as CP
import math


def finding_props(temp, fluid):
    fluid_Props = dict()

    fluid_Props['density'] = CP.PropsSI('D', 'T', temp, 'P', 101325, fluid)

    fluid_Props['Pr'] = CP.PropsSI(
        'PRANDTL', 'T', temp, 'P', 101325, fluid)

    fluid_Props['Mu'] = CP.PropsSI(
        'VISCOSITY', 'T', temp, 'P', 101325, fluid)

    fluid_Props['Cp'] = CP.PropsSI('CPMASS', 'T', temp, 'P', 101325, fluid)

    fluid_Props['k'] = CP.PropsSI('L', 'T', temp, 'P', 101325, fluid)

    return fluid_Props


def calculating_h(m, di, do, Di, fluid_Props, fluid):
    Rho = fluid_Props['density']
    k = fluid_Props['k']
    Mu = fluid_Props['Mu']
    Pr = fluid_Props['Pr']

    # hot fluid
    if fluid == 'hydrogen':
        # calculating Vc
        Vc = m / (Rho * (math.pi/4) * di**2)
        # claculating Re
        Re = (Rho * Vc * di) / Mu
    # cold fluid
    if fluid == 'ethanol':
        # calculating Vc
        Vc = m / (Rho * (math.pi/4) * (Di**2 - do**2))
        # claculating Re
        Re = (Rho * Vc * (Di-do)) / Mu

    # calculating Nu
    if Re > 2300:
        if fluid == 'hydrogen':
            # cooling of fluid
            n = 0.3
            Nu = 0.023*(Re**0.8)*(Pr**n)
        if fluid == 'ethanol':
            # heating of fluid
            n = 0.4
            Nu = 0.023*(Re**0.8)*(Pr**n)
    if Re < 2300:
        Nu = 4.36
        # Nu = 3.66

    # calculating h
    # hot fluid
    if fluid == 'hydrogen':
        h = Nu * k / di
    # cold fluid
    if fluid == 'ethanol':
        d = (Di**2 - do**2)/do
        h = Nu * k / d

    return h


def calculating_Q(m_h, m_c, props_hot, props_cold, E, Th, Tc):
    Cp_hot = props_hot['Cp']
    Cp_cold = props_cold['Cp']

    mCp_min = min(m_h*Cp_hot, m_c*Cp_cold)
    Q_max = mCp_min*(Th-Tc)
    Q = E*Q_max
    return Q


def calculating_outlet_temp(Q, m, inlet_temp, fluid_props, fluid):
    C_p = fluid_props['Cp']

    if fluid == 'hydrogen':
        outlet_temp = inlet_temp - Q/(m*C_p)
    if fluid == 'ethanol':
        outlet_temp = inlet_temp + Q/(m*C_p)

    return outlet_temp


def calculating_effectiveness(Th_i, Th_o, Tc_i, Tc_o, m_c, m_h, h_hot, h_cold, di, do, Q, props_hot, props_cold):
    # Constant values
    k_pipe = 15.1
    Rfh = 8.815E-05
    Rfc = 0.0003526
    Cp_hot = props_hot['Cp']
    Cp_cold = props_cold['Cp']

    # calculating U
    U = (
        1/h_cold + Rfc + (di/2*k_pipe)*math.log(do/di) + (do/di)*Rfh + (do/di)*(1/h_hot))**(-1)

    # calculating LMTD
    LMTD = ((Th_i-Tc_o)-(Th_o-Tc_i))/math.log((Th_i-Tc_o)/(Th_o-Tc_i))

    # calculating Area
    area = Q/(U*LMTD)

    # calculating NTU
    NTU = (U*area)/min(m_h*Cp_hot, m_c*Cp_cold)

    # calculating C
    C = min(m_h*Cp_hot, m_c*Cp_cold)/max(m_h*Cp_hot, m_c*Cp_cold)

    # calculating effectiveness
    E = (1-math.exp(NTU*(C-1)))/(1-C*math.exp(NTU*(C-1)))

    print("Th_o:", Th_o-273.15)
    print("Tc_o:", Tc_o-273.15)
    print("LMTD:", LMTD)
    print("Area:", area)
    print("NTU:", NTU)

    return E


def HX_Designing(di, do, Di, Th, Tc, m_h, m_c, E, hot_fluid, cold_fluid):
    props_hot = finding_props(Th, hot_fluid)
    props_cold = finding_props(Tc, cold_fluid)

    h_hot = calculating_h(m_h, di, do, Di, props_hot, hot_fluid)
    h_cold = calculating_h(m_c, di, do, Di, props_cold, cold_fluid)

    Q = calculating_Q(
        m_h, m_c, props_hot, props_cold, E, Th, Tc)

    Th_o = calculating_outlet_temp(Q, m_h, Th, props_hot, hot_fluid)
    Tc_o = calculating_outlet_temp(Q, m_c, Tc, props_cold, cold_fluid)

    new_E = calculating_effectiveness(
        Th, Th_o, Tc, Tc_o, m_c, m_h, h_hot, h_cold, di, do, Q, props_hot, props_cold)

    if (Tc_o-Tc <= 0.01):
        print("Heat Tranfer:", Q)
        return new_E
    else:
        print("Heat Tranfer:", Q)
        print("Heat Exchanger Effectiveness:", new_E)
        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
        print("\n")
        return HX_Designing(
            di, do, Di, (Th+Th_o)/2, (Tc+Tc_o)/2, m_h, m_c, E, hot_fluid, cold_fluid)


#
#
#
#
#
# Inputs
di = 3.388e-03
do = 3.81e-03
Di = 6.35e-03
Th_i = 150 + 273.15
Tc_i = 25 + 273.15
m_h = 0.5
m_c = 0.3

hot_fluid = 'hydrogen'
cold_fluid = 'ethanol'

assumed_E = 0.75

calc_E = HX_Designing(
    di, do, Di, Th_i, Tc_i, m_h, m_c, assumed_E, hot_fluid, cold_fluid)

print("Final Heat Exchanger Effectiveness:", calc_E)
