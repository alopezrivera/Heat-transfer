import time
import datetime
import numpy as np

from src import constants as c


class Tank():
    def __init__(self):
        self.k_c = c.k_cf
        self.t_c = c.t_cf
        self.k_a = c.k_alu
        self.t_a = c.t_alu
        self.A   = c.S
        self.m   = c.m_tot  # kg
        self.V   = c.V_tot  # m3


class Lamp():
    def __init__(self):
        self.eps          = c.emissivity
        self.T            = c.lamp_temp
        self.shape_factor = c.shape_factor


# dummy class
class N2O():
    def __init__(self):
        self.rho = c.N20_rho
        self.cp  = c.N20_cp


class Thermals():
    def __init__(self, tank, lamp, n2o, T, dt, T_amb, evap=False):
        self.tank  = tank
        self.lamp  = lamp
        self.n2o   = n2o
        self.T     = T
        self.dt    = dt
        self.T_amb = T_amb
        self.evap  = evap

    def get_fluid_properties(self, T):
        vap_frac = self.n2o.vap_mass_frac(T)  # mass fraction of vapor for mass averaged density and Cp
        # mass_frac = self.n2o.m_V(T)/self.n2o.m_L(T)
        # print(vap_frac)
        rho = (1 - vap_frac) * self.n2o.rho_V(T) + vap_frac * self.n2o.rho_L(T)
        cp  = (1 - vap_frac) * self.n2o.cp_V(T) + vap_frac * self.n2o.cp_L(T)

        return rho, cp, vap_frac

    def P_in(self, T, T_lamp):
        h_r = 4 * self.lamp.eps * c.sigma * self.lamp.shape_factor * ((T + self.lamp.T) / 2) ** 3
        R = (1 / h_r + self.tank.t_c / self.tank.k_c + self.tank.t_a / self.tank.k_a) * 1 / self.tank.A

        dT = T_lamp - T

        return dT/R

    def Q_radiated(self):
        pass

    def forward_euler(self, cutoff):

        t0 = time.time()

        self.time     = np.arange(0, self.T, self.dt)
        self.T_tank   = np.ones(len(self.time)) * self.T_amb
        self.rho_n2o  = np.ndarray(len(self.time))
        self.cp_n2o   = np.ndarray(len(self.time))
        self.vap_frac = np.ndarray(len(self.time))

        self.T_tank[0] = self.n2o.T0
        self.rho_n2o[0], self.cp_n2o[0], self.vap_frac[0] = self.get_fluid_properties(self.T_tank[0])

        def f(T_tank, vap_frac0,
              T_lamp, tank, n2o,
              P_in,
              dt,
              evap,
              log=False
              ):

            rho_n2o1, cp_n2o1, vap_frac1 = self.get_fluid_properties(T_tank)

            # Incoming power
            p_in = P_in(T=T_tank, T_lamp=T_lamp)

            # Radiated power
            p_rad = dt*c.sigma*((T_tank)**4-c.T_amb**4)*c.S

            # Vaporization power
            evap_mass = (vap_frac1 - vap_frac0) * tank.m
            p_evap = dt*n2o.hvap * evap_mass if evap else 0

            dT_tank = (p_in - p_evap - p_rad) / (rho_n2o1 * cp_n2o1 * tank.V)

            if log:
                print(f'{(p_in - p_evap - p_rad):.3f} [kW]')
                print(f'{dT_tank:.3f} [K]')

            return dT_tank

        for i in range(1, len(self.time)):

            print(f'%{(i/len(self.time)*100):.2f}')

            self.T_tank[i] = self.T_tank[i-1] + self.dt*f(T_tank=self.T_tank[i-1],
                                                          vap_frac0=self.vap_frac[i - 1],
                                                          T_lamp=self.lamp.T,
                                                          P_in=self.P_in,
                                                          dt=self.dt,
                                                          tank=self.tank,
                                                          n2o=self.n2o,
                                                          evap=self.evap)

            if self.n2o.p(self.T_tank[i]) > cutoff:

                runtime = time.time()

                print('max pressure reached after  :: ', f'[h:m:s]  {str(datetime.timedelta(seconds=self.time[i]))}')
                print('run time                    :: ', f'[s]      {runtime - t0:.2f}')

                break
