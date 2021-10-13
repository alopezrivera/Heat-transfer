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
        self.T            = c.T_lamp
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
        vap_frac = self.n2o.vap_mass_frac(T)    # mass fraction of vapor for mass averaged density and Cp
        rho = (1 - vap_frac) * self.n2o.rho_V(T) + vap_frac * self.n2o.rho_L(T)
        cp  = (1 - vap_frac) * self.n2o.cp_V(T) + vap_frac * self.n2o.cp_L(T)

        return rho, cp, vap_frac

    def R_rad(self, T_cf0, f_12=c.epsilon):
        """
        Heat transfer coefficient and thermal resistance:

            h_r = sigma * f * (T_1 - T_2)**4/(T_1 - T_2)

            R   = 1/(h_r*A_1)

        where
            - sigma: Boltzmann constant
            - f    : Reflection transfer function  ->  epsilon   if 1 (reflecting object) is surrounded by 2
                                                   ->  1         if 1 is a black body
            - A_1  : Area of radiating (emitting + reflecting) body
            - T_1  : Surface temperature of radiating body
            - T_2  : Temperature of surroundings

        Assumptions:
            1. The tank radiates heat equally over its entire surface
                - Clearly wrong
                    - Reflection is clearly higher directly under the lamps
                    - The temperature distribution of the tank is not homogenous -> emission cannot be
                - HOWEVER
                    - The heat transfer through the Alu and CF is considered to happen through the area illuminated
                      by the lamps, and temperature of the tank is considered homogenous
                        - The temperature of the tank will that of the illuminated patches -> radiation will be high
                        - The assumption is considered valid as heat is transferred to the gas through the illuminated
                          spots
                            - Tradeoff is between realistic heat transfer to the gas and realistic tank temperature
                                - Decision in favor of realistic gas heat influx
                                - Overestimated tank temperature will cause higher radiation = heat rejection
                                    - Overestimation of the tank's heating time, which is considered acceptable
        """
        if T_cf0 - c.T_amb > 0:
            h_r = c.sigma * f_12 * (T_cf0 - c.T_amb)**4/(T_cf0 - c.T_amb)
            return 1/(h_r*c.S)
        else:
            return 0

    def R_cf(self):
        """
        Assumptions:
            - Heat transfer happens through the area illuminated by the lamps
                - Thermal resistance will be much lower
        """
        L = c.t_cf
        k = c.k_cf
        A = c.s_lamp
        return L/(k*A)

    def R_alu(self):
        """
        Assumptions:
            - Heat transfer happens through the area illuminated by the lamps
                - Thermal resistance will be much lower
        """
        L = c.t_alu
        k = c.k_alu
        A = c.s_lamp
        return L/(k*A)

    def R_tot(self, T_cf0):
        """
        Total thermal resistance
        """
        return self.R_rad(T_cf0) + self.R_cf() + self.R_alu()

    def Q_dot(self, T_alu1, T_cf0):
        """
        Heat flux
        """
        return (self.lamp.T - T_alu1)/self.R_tot(T_cf0)

    def dT(self,
           T_tank,
           Q_dot,
           vap_frac0,
           log):
        """
        Temperature raise

        :param T_tank:      Tank temperature
        :param vap_frac0:   Nitrous vapour fraction
        :param log:         Log heat flux and temperature raise
        """
        # Fluid properties
        rho_n2o1, cp_n2o1, vap_frac1 = self.get_fluid_properties(T_tank)

        # Vaporization power
        if self.evap:
            evap_mass = (vap_frac1 - vap_frac0) * self.tank.m
            Q_dot_evap = self.dt * self.n2o.hvap * evap_mass
        else:
            Q_dot_evap = 0

        # Temperature raise
        dT_tank = (Q_dot - Q_dot_evap) / (rho_n2o1 * cp_n2o1 * self.tank.V)

        if log:
            print(f'{(Q_dot - Q_dot_evap):.3f} [kW]')
            print(f'{dT_tank:.3f} [K]')

        return dT_tank

    def simulate(self, cutoff):
        """
        Conduct heating simulation
        """
        self.t        = np.arange(0, self.T, self.dt)

        # Temperatures
        self.T_cf0    = self.T_amb  * np.ones(len(self.t))
        self.T_alu0   = self.n2o.T0 * np.ones(len(self.t))
        self.T_alu1   = self.n2o.T0 * np.ones(len(self.t))
        self.T_n2o    = self.n2o.T0 * np.ones(len(self.t))

        # Gas characteristics
        self.rho_n2o, self.cp_n2o, self.vap_frac = [prop * np.ones(len(self.t)) for prop in self.get_fluid_properties(self.n2o.T0)]

        # Time
        t0 = time.time()

        for i in range(len(self.t)):

            print(f'%{(i / len(self.t) * 100):.2f}')

            Q_dot = self.Q_dot(T_cf0=self.T_cf0[i-1], T_alu1=self.T_alu1[i-1])

            # Resistances
            R_rad = self.R_rad(self.T_cf0[i-1])
            R_cf  = self.R_cf()
            R_alu = self.R_alu()

            # Carbon fiber temperature
            self.T_cf0[i]  = self.lamp.T - Q_dot*(R_rad)
            self.T_alu0[i] = self.lamp.T - Q_dot*(R_rad + R_cf)
            self.T_alu1[i] = self.lamp.T - Q_dot*(R_rad + R_cf + R_alu)

            # Forward Euler integration scheme
            self.T_n2o[i] = self.T_n2o[i - 1] + self.dT(Q_dot=Q_dot,
                                                        T_tank=self.T_n2o[i-1],
                                                        vap_frac0=self.vap_frac[i-1],
                                                        log=False)

            if self.n2o.p(self.T_n2o[i]) > cutoff:

                runtime = time.time()

                print('max pressure reached after  :: ', f'[h:m:s]  {str(datetime.timedelta(seconds=self.t[i]))}')
                print('run time                    :: ', f'[s]      {runtime - t0:.2f}')

                break

    # """
    # Previous model
    # """
    # def P_in(self, T, T_lamp):
    #     h_r = 4 * self.lamp.eps * c.sigma * self.lamp.shape_factor * ((T + self.lamp.T) / 2) ** 3
    #     R = (1 / h_r + self.tank.t_c / self.tank.k_c + self.tank.t_a / self.tank.k_a) * 1 / self.tank.A
    #
    #     dT = T_lamp - T
    #
    #     return dT/R
    #
    # def forward_euler(self, cutoff):
    #
    #     t0 = time.time()
    #
    #     self.t     = np.arange(0, self.T, self.dt)
    #     self.T_n2o   = np.ones(len(self.t)) * self.T_amb
    #     self.rho_n2o  = np.ndarray(len(self.t))
    #     self.cp_n2o   = np.ndarray(len(self.t))
    #     self.vap_frac = np.ndarray(len(self.t))
    #
    #     self.T_n2o[0] = self.n2o.T0
    #     self.rho_n2o[0], self.cp_n2o[0], self.vap_frac[0] = self.get_fluid_properties(self.T_n2o[0])
    #
    #     def f(T_tank, vap_frac0,
    #           T_lamp, tank, n2o,
    #           P_in,
    #           dt,
    #           evap,
    #           log=False
    #           ):
    #
    #         # Incoming power
    #         p_in = P_in(T=T_tank, T_lamp=T_lamp)
    #
    #         # Radiated power
    #         p_rad = dt*c.sigma*((T_tank)**4-c.T_amb**4)*c.S
    #
    #         # Vaporization power
    #         rho_n2o1, cp_n2o1, vap_frac1 = self.get_fluid_properties(T_tank)
    #         evap_mass = (vap_frac1 - vap_frac0) * tank.m
    #         p_evap = dt*n2o.hvap * evap_mass if evap else 0
    #
    #         dT_tank = (p_in - p_evap - p_rad) / (rho_n2o1 * cp_n2o1 * tank.V)
    #
    #         if log:
    #             print(f'{(p_in - p_evap - p_rad):.3f} [kW]')
    #             print(f'{dT_tank:.3f} [K]')
    #
    #         return dT_tank
    #
    #     for i in range(1, len(self.t)):
    #
    #         print(f'%{(i/len(self.t)*100):.2f}')
    #
    #         self.T_n2o[i] = self.T_n2o[i-1] + self.dt*f(T_tank=self.T_n2o[i-1],
    #                                                       vap_frac0=self.vap_frac[i - 1],
    #                                                       T_lamp=self.lamp.T,
    #                                                       P_in=self.P_in,
    #                                                       dt=self.dt,
    #                                                       tank=self.tank,
    #                                                       n2o=self.n2o,
    #                                                       evap=self.evap)
    #
    #         if self.n2o.p(self.T_n2o[i]) > cutoff:
    #
    #             runtime = time.time()
    #
    #             print('max pressure reached after  :: ', f'[h:m:s]  {str(datetime.timedelta(seconds=self.t[i]))}')
    #             print('run time                    :: ', f'[s]      {runtime - t0:.2f}')
    #
    #             break
