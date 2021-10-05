import numpy as np

# Time
dt = 0.00001                    # s
t = 0                           # s

# Temperature
T_amb = 20.                     # C
T_amb = T_amb+273.15
T_end = 301.                    # K

# Heat
eff_HL = 0.7                    # -
n_HL = 8                        # -
P_HL = 650.                     # W
Q_in = eff_HL*n_HL*P_HL*dt      # J

# Dimensions
L = 4.7                         # m
d = 0.2778                      # m
S = np.pi*d*L                   # m2

# Mass and volume
m_tot = 174.                    # kg
V_tot = 0.259                   # m3

# Aluminium
m_al = 16.7                     # kg
k_al = 237.                     # W/m K
cp_al= 890.                     # W/kgK
t_al = 0.0015                   # m
T_al = [T_amb]                  # K
q_al = [0.]                     # W/m2

# Carbon fiber
m_cf = 12.9                     # kg
k_cf = 6.                       # W/m K 5-7
cp_cf= 800.                     # W/kgK
t_cf = 0.0025                   # m
T_cf = [T_amb]                  # K
q_cf = [0.]                     # W/m2
