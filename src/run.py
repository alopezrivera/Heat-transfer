import numpy as np
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
emissivity = 0.01
lamp_temp = 350         # K

# sim data 
T = 1000                # s
dt = 0.01             # s     0.001 for stability if evap is considered
T0 = 288                # K


tank = Tank(k_carbon, t_carbon, k_alu, t_alu, A_tank)
lamp = Lamp(emissivity, lamp_temp)
n2o = N2O()             


thermals = Thermals(tank, lamp, n2o, T, dt, T0)
thermals.runge_kutta4(p_max)

plt.plot(thermals.time, thermals.T_tank)
plt.xlabel('time [s]')
plt.ylabel('temp [K]')
plt.show()
