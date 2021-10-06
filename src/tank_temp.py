import numpy as np
from matplotlib import pyplot as plt


class Tank():
    def __init__(self, k_carbon, t_carbon, k_alu, t_alu, A):
        self.k_c = k_carbon
        self.t_c = t_carbon
        self.k_a = k_alu
        self.t_a = t_alu
        self.A = A
        self.m = 174       # kg
        self.V = 0.259      # m3


class Lamp():
    def __init__(self, eps, Temp):
        self.eps = eps
        self.T = Temp
        self.sig = 5.67e-8


# dummy class 
class N2O():
    def __init__(self):
        self.rho = 300
        self.cp = 1000


class Thermals():
    def __init__(self, tank, lamp, n2o, T, dt, T_amb):
        self.tank = tank
        self.lamp = lamp
        self.n2o = n2o
        self.T = T
        self.dt = dt 
        self.T_amb = T_amb

    def thermal_res(self, T):
        h_r = 4 * self.lamp.eps * self.lamp.sig * ((T + self.lamp.T)/2)**3
        R = (1/h_r + self.tank.t_c/self.tank.k_c + self.tank.t_a/self.tank.k_a) * 1/self.tank.A
        
        return R

    def get_fluid_properties(self, T):
        vap_frac = self.n2o.vap_mass_frac(T)            # mass fraction of vapor
        #mass_frac = self.n2o.m_V(T)/self.n2o.m_L(T)
        #print(vap_frac, mass_frac)
        rho = (1-vap_frac) * self.n2o.rho_V(T) + vap_frac * self.n2o.rho_L(T) 
        cp = (1-vap_frac) * self.n2o.cp_V(T) + vap_frac * self.n2o.cp_L(T) 
        
        return rho, cp, vap_frac

    def runge_kutta4(self, cutoff):
        self.time = np.arange(0, self.T, self.dt)

        self.T_tank = np.ones(len(self.time))
        self.dT_tank = np.ones(len(self.time))

        self.rho_n2o = np.ndarray(len(self.time))
        self.cp_n2o = np.ndarray(len(self.time))
        self.vap_frac = np.ndarray(len(self.time))

        # Initialization
        self.T_tank *= self.T_amb
        self.rho_n2o[0], self.cp_n2o[0], self.vap_frac[0] = self.get_fluid_properties(self.T_tank[0])

        def f(T_tank, vap_frac0,
              thermal_res,
              T_lamp, tank, n2o
              ):

            R         = thermal_res(T_tank)
            deltaT    = T_lamp - T_tank

            rho_n2o1, cp_n2o1, vap_frac1 = self.get_fluid_properties(T_tank)

            evap_mass = (vap_frac1 - vap_frac0) * tank.m
            Q_evap    = n2o.hvap * evap_mass

            dT_tank   = (deltaT/R - Q_evap) / (rho_n2o1*cp_n2o1*tank.V)

            return dT_tank

        for i in range(1, len(self.time)):

            print(f'{(i / len(self.time) * 100):.2f}%')

            k1a = self.dt * self.dT_tank[i-1]
            k1b = self.dt * f(T_tank=self.T_tank[i - 1],
                              vap_frac0=self.vap_frac[i - 1],
                              T_lamp=self.lamp.T,
                              thermal_res=self.thermal_res,
                              tank=self.tank,
                              n2o=self.n2o)

            k2a = self.dt * (self.dT_tank[i-1]+k1b/2)
            k2b = self.dt * f(T_tank=self.T_tank[i - 1] + k1a/2,
                              vap_frac0=self.vap_frac[i - 1],
                              T_lamp=self.lamp.T,
                              thermal_res=self.thermal_res,
                              tank=self.tank,
                              n2o=self.n2o)

            k3a = self.dt * (self.dT_tank[i-1]+k2b/2)
            k3b = self.dt * f(T_tank=self.T_tank[i - 1] + k2a/2,
                              vap_frac0=self.vap_frac[i - 1],
                              T_lamp=self.lamp.T,
                              thermal_res=self.thermal_res,
                              tank=self.tank,
                              n2o=self.n2o)

            k4a = self.dt * (self.dT_tank[i-1]+k3b)
            k4b = self.dt * f(T_tank=self.T_tank[i - 1] + k3a,
                              vap_frac0=self.vap_frac[i - 1],
                              T_lamp=self.lamp.T,
                              thermal_res=self.thermal_res,
                              tank=self.tank,
                              n2o=self.n2o)

            self.T_tank[i]  = self.T_tank[i-1]  + 1/6*(k1a + k2a/2 + k3a/2 + k4a)
            self.dT_tank[i] = self.dT_tank[i-1] + 1/6*(k1b + k2b/2 + k3b/2 + k4b)

            print('================')
            print('Euler')
            print(self.dt * f(T_tank=self.T_tank[i - 1],
                              vap_frac0=self.vap_frac[i - 1],
                              T_lamp=self.lamp.T,
                              thermal_res=self.thermal_res,
                              tank=self.tank,
                              n2o=self.n2o))
            print('================')
            print('Runge Kutta')
            print(1/6*(k1a + k2a/2 + k3a/2 + k4a))
            print('================')
            # print(self.dT_tank[i])
            print('\n\n')

            self.rho_n2o[i], self.cp_n2o[i], self.vap_frac[i] = self.get_fluid_properties(self.T_tank[i-1])

            if self.n2o.p(self.T_tank[i]) > cutoff:
                print('max pressure reached after: ', self.time[i], ' seconds')
                break

    def forward_euler(self, cutoff):
        self.time = np.arange(0, self.T, self.dt)
        self.T_tank = np.ones(len(self.time))*288
        self.rho_n2o = np.ndarray(len(self.time))
        self.cp_n2o = np.ndarray(len(self.time))
        self.vap_frac = np.ndarray(len(self.time))

        self.T_tank[0] = self.T_amb
        self.rho_n2o[0], self.cp_n2o[0], self.vap_frac[0] = self.get_fluid_properties(self.T_tank[0])

        for i in range(1, len(self.time)):

            print(f'{(i/len(self.time)*100):.2f}%')

            R = self.thermal_res(self.T_tank[i-1])
            deltaT = self.lamp.T - self.T_tank[i-1]

            self.rho_n2o[i], self.cp_n2o[i], self.vap_frac[i] = self.get_fluid_properties(self.T_tank[i-1])
            
            evap_mass = (self.vap_frac[i] - self.vap_frac[i-1]) * self.tank.m
            Q_evap = self.n2o.hvap * evap_mass

            self.T_tank[i] = self.T_tank[i-1] + self.dt * ((deltaT/R - Q_evap) / self.rho_n2o[i] / self.cp_n2o[i] / self.tank.V)

            if self.n2o.p(self.T_tank[i]) > cutoff:
                print('max pressure reached after: ', self.time[i], ' seconds')
                break


if __name__ == '__main__':

    tank = Tank(20, 2e-3, 200, 3e-3, 4)
    lamp = Lamp(0.6, 400)
    n2o = N2O()
    thermals = Thermals(tank, lamp, n2o, 3600, 10, 288)
    thermals.forward_euler()

    plt.plot(thermals.time, thermals.T_tank)
    plt.xlabel('time [s]')
    plt.ylabel('temp [K]')
    plt.show()
