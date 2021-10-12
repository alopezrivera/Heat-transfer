from mpl_plotter.two_d import line

from src.tank_temp import Tank, Lamp, Thermals
from src import constants as c
from src.n2o_class import N2O

from alexandria.math.differentiation import derivative


# Setup
T     = 2*60**2       # s
dt    = 0.1           # s  0.001 for stability if evap is considered
T_amb = 288           # K


if __name__ == '__main__':

    tank = Tank()
    lamp = Lamp()
    n2o = N2O()

    thermals = Thermals(tank=Tank(),
                        lamp=Lamp(),
                        n2o=N2O(),
                        T=T,
                        dt=dt,
                        T_amb=T_amb,
                        evap=False)

    thermals.forward_euler(c.p_max)

    line(thermals.time/60,
         thermals.T_tank,
         # thermals.n2o.p(thermals.T_tank)/100,
         # derivative(thermals.time/60/60, thermals.n2o.p(thermals.T_tank)/100),
         # x_label='Time [min]',
         # y_label='Temp [K]',
         # y_bounds=[0, 60],
         y_bounds=[280, 302],
         color='C2',
         show=True
         )
