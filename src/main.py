# Written in an afternoon on launch campaign, please check and use at own risk.
# The programming style isn't pretty nor fast but seems to do its job.
# Rolf Wubben 03/10/2021

import numpy as np
from matplotlib import pyplot as plt

from src.methods import p, rho_V, rho_L, cp_V, cp_L, h_V, h_L, Q_rad

from src.constants import m_tot, V_tot, T_end, t, T_cf, Q_in, cp_cf, m_cf, S, dt


T_N2O = [271.53677] #K
p_N2O = [p(T_N2O[-1])] #kPa
rho_V_N2O = [rho_V(T_N2O[-1])]
rho_L_N2O = [rho_L(T_N2O[-1])]
cp_V_N2O = [cp_V(T_N2O[-1])]
cp_L_N2O = [cp_L(T_N2O[-1])]
h_V_N2O = [h_V(T_N2O[-1])]
h_L_N2O = [h_L(T_N2O[-1])]
a_N2O = [(m_tot/V_tot-rho_V_N2O[-1])/(rho_L_N2O[-1]-rho_V_N2O[-1])]
m_V_N2O = [(1-a_N2O[-1])*rho_V_N2O[-1]*V_tot]
m_L_N2O = [a_N2O[-1]*rho_L_N2O[-1]*V_tot]
Q_vap = [0.]

while T_N2O[-1] < T_end and t < 10000:
    a_N2O   = (m_tot / V_tot - rho_V_N2O[-1]) / (rho_L_N2O[-1] - rho_V_N2O[-1])
    m_V_N2O = (1 - a_N2O[-1]) * rho_V_N2O[-1] * V_tot
    m_L_N2O = a_N2O[-1] * rho_L_N2O[-1] * V_tot

    dm = (m_L_N2O[-2]-m_L_N2O[-1])
    dh = h_V_N2O[-1]-h_L_N2O[-1]

    dH = dm*dh

    dp = (p_N2O[-2]-p_N2O[-1])

    Q_vap   = dH - dp/rho_V_N2O[-2]*dm


# Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     print_hi('PyCharm')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/

