import numpy as np
from pynverse import inversefunc

# Temperature
T_amb = 20+273.15               # K
T_end = 300.8                   # K

# Dimensions
L = 4.7                         # m
d = 0.2778                      # m
S = np.pi*d*L                   # m2

# Maximum pressure
p_max = 6000                    # kpa

# Mass and volume
m_tot = 174.                    # kg
V_tot = 0.259                   # m3

# Nitrous
T_c = 309.57                    # K
p_c = 7251.                     # kPa

B = [-6.71893,
     1.35966,
     -1.3779,
     -4.051]
N = [1,
     3/2,
     5/2,
     5]

p = lambda T: p_c * np.exp(1
                           / (T/T_c)
                           * sum(
                                    [
                                        B[i] * (1-T/T_c) ** N[i] for i in range(len(B))
                                    ]
                                )
                           )

N2O_T0 = 280.425                # K

# Aluminium
m_alu  = 16.7                   # kg
k_alu  = 237.                   # W/m/K
cp_alu = 890.                   # W/kg/K
t_alu  = 0.0015                 # m
T_alu  = [N2O_T0]               # K
q_alu  = [0.]                   # W/m2

# Carbon fiber
m_cf  = 12.9                    # kg
k_cf  = 6.                      # W/m K 5-7
cp_cf = 800.                    # W/kgK
t_cf  = 0.0025                  # m
T_cf  = [T_amb]                 # K
q_cf  = [0.]                    # W/m2

# Lamp data

n_lamps      = 2

emissivity   = 0.6
T_lamp       = 350                              # K
L_lamp       = 0.73                             # m     Length of the lamp
W_lamp       = 0.05                             # m     "Width" of the lamp
D_lamp       = 0.1                              # m     Distance from cylinder
c_lamp       = D_lamp                           # m     Area illuminated by the lamp: cord of the tank
s_lamp       = 2*d/2*np.arcsin(c_lamp/2/d/2)    # m     Area illuminated by the lamp: distance along surface of cylinder
S_lamp       = L_lamp*s_lamp*n_lamps            # m^2   Area illuminated by the lamp: area

shape_factor = (c_lamp/2/d/2)**2                # -     sin(theta)**2 [Mills and Coimbra, Table 6.1, entry 12]

# Radiation
sigma        = 5.67e-8          # W/m^2/K^4
epsilon      = 0.5
