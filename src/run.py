import time
from copy import deepcopy
import matplotlib as mpl
mpl.use('Qt5Agg')
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
T = 100                # s
dt = 0.001              # s  0.001 for stability if evap is considered
T0 = 288                # K


if __name__ == '__main__':

    tank = Tank(k_carbon, t_carbon, k_alu, t_alu, A_tank)
    lamp = Lamp(emissivity, lamp_temp, shape_factor)
    n2o = N2O()

    t0 = time.time()
    thermals = Thermals(tank, lamp, n2o, T, dt, T0, evap=True)

    t1, t2 = deepcopy(thermals), deepcopy(thermals)

    t1.forward_euler(p_max)
    runtime1 = time.time()

    print('run time: ', runtime1 - t0, ' seconds')

    t2.runge_kutta4(p_max)
    runtime2 = time.time()

    print('run time: ', runtime2 - t0, ' seconds')

    plt.plot(t1.time, t1.T_tank)
    plt.plot(t2.time, t2.T_tank)
    plt.xlabel('time [s]')
    plt.ylabel('temp [K]')
    plt.show()
