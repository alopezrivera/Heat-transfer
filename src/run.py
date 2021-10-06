import numpy as np
from copy import deepcopy
from matplotlib import pyplot as plt

from src.tank_temp import Tank, Lamp, Thermals
from src.n2o_class import N2O

# tank data
k_carbon = 10           # W/m/K
k_alu = 200             # W/m/K
t_carbon = 2e-3         # m
t_alu = 3e-3            # m
A_tank = 4              # m^2
p_max = 7000            # kpa

# lamp data
shape_factor = 0.00269
emissivity = 0.6
lamp_temp = 350         # K

# sim data
T = 3600                # s
dt = 0.001              # s  0.001 for stability if evap is considered
T0 = 288                # K


tank = Tank(k_carbon, t_carbon, k_alu, t_alu, A_tank)
lamp = Lamp(emissivity, lamp_temp)
n2o = N2O()             


thermals = Thermals(tank, lamp, n2o, T, dt, T0)

t1, t2 = deepcopy(thermals), deepcopy(thermals)
t1.forward_euler(p_max)
t2.runge_kutta4(p_max)

plt.plot(t1.time, t1.T_tank)
plt.plot(t2.time, t2.T_tank)
plt.xlabel('time [s]')
plt.ylabel('temp [K]')
plt.show()
