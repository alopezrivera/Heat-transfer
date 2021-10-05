import numpy as np

from src.constants import S, T_amb, dt


def p(T):
    T_c = 309.57             # K
    p_c = 7251.              # kPa
    T_r = T / T_c
    a = 1 - T_r
    b = [-6.71893, 1.35966, -1.3779, -4.051]
    n = [1., 3./2., 5./2., 5.]
    p = 0
    for i in range(len(b)):
        p += b[i]*a**n[i]
    return p_c*np.exp(1/T_r*p)


def rho_L(T):
    T_c = 309.57             # K
    rho_c = 452.             # kg/m3
    T_r = T / T_c
    a = 1 - T_r
    b = [1.72328, -0.83950, 0.51060, -0.10412]
    n = [1./3., 2./3., 3./3., 4./3.]
    rho = 0
    for i in range(len(b)):
        rho += b[i]*a**n[i]
    return rho_c*np.exp(rho)


def rho_V(T):
    T_c = 309.57             # K
    rho_c = 452.             # kg/m3
    T_r = T / T_c
    a = 1/T_r-1
    b = [-1.00900, -6.28792, 7.50332, -7.90463, 0.629427]
    n = [1./3., 2./3., 3./3., 4./3., 5/3]
    rho = 0
    for i in range(len(b)):
        rho += b[i]*a**n[i]
    return rho_c*np.exp(rho)


def h_L(T):
    T_c = 309.57             # K
    T_r = T / T_c
    a = 1 - T_r
    b = [-200., 116.043, -917.225, 794.779, -589.587]
    n = [0, 1./3., 2./3., 3./3., 4./3.]
    h = 0
    for i in range(len(b)):
        h += b[i]*a**n[i]
    return h


def h_V(T):
    T_c = 309.57             # K
    T_r = T / T_c
    a = 1 - T_r
    b = [-200., 440.055, -459.701, 434.081, -485.338]
    n = [0, 1./3., 2./3., 3./3., 4./3.]
    h = 0
    for i in range(len(b)):
        h += b[i]*a**n[i]
    return h


def cp_L(T):
    T_c = 309.57             # K
    T_r = T / T_c
    a = 1 - T_r
    b = [2.49973, 0.023454, -3.80136, 13.0945, -14.5180]
    n = [0, -1., 1, 2, 3]
    cp = 1
    for i in range(1,len(b)):
        cp += b[i]*a**n[i]
    return b[0]*cp


def cp_V(T):
    T_c = 309.57             # K
    T_r = T / T_c
    a = 1 - T_r
    b = [132.632, 0.052187, -0.364923, -1.20233, 0.536141]
    n = [0, -2./3., -1./3., 1./3., 2./3.]
    cp = 1
    for i in range(1,len(b)):
        cp += b[i]*a**n[i]
    return b[0]*cp


def Q_rad(T):
    # S = np.pi*0.278*4.7
    sb = 5.67e-8            # Wm-2K-4
    return sb*S*(T**4-T_amb**4)*dt