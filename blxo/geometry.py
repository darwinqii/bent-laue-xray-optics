import os
file_dir = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.extend([file_dir])
import numpy as np
from scipy.optimize import fsolve


class Angles:
    def __init__(self, chi, theta, nu, t, r, p):
        self.chi = chi
        self.theta = theta
        self.nu = nu
        self.t = t
        self.r = r
        self.p = p

    # magic condition core angles
    def delta_psi_AB(self):
        return -self.t * np.tan(self.theta - self.chi) / self.r

    def delta_chi_AB(self):
        return -(1 + self.nu) * self.t * np.sin(2 * self.chi) / (2 * self.r)

    def alpha_AB(self):
        return (self.t / self.p) * (np.sin(2 * self.theta) / np.cos(self.chi - self.theta))

    def delta_theta_AB(self):
        '''
        Notes:
            This $\Delta\theta_{AB}$ is caused by the finite source distance.
            ONLY when the 'magic condition' is met, this is the difference of Bragg angles between point A and B.
            When the 'magic condition' is not met, a misalignment angle should be added for the difference of Bragg
            angles between A and B.
        Returns:

        '''
        return -1 / 2 * self.alpha_AB()

    # Other alpha angles (angles from the finite source distance)
    def alpha_BC(self):  # todo: double check using point A or point C for theta and chi
        qmb_foot_length = Lengths(self.chi, self.theta, self.nu, self.t, self.r, self.p).foot_length()
        return np.cos(self.theta_B() + self.chi_B()) * qmb_foot_length / self.p

    def chi_B(self):
        return self.chi + self.delta_chi_AB()

    def chi_C(self):
        return self.chi_B()

    def theta_misalign(self):
        try:
            from blxo.mc import magic_condition_angles
        except:
            from mc import magic_condition_angles
        return magic_condition_angles(self.chi, self.theta, self.nu, self.t, self.r, self.p)

    def theta_B(self):
        return self.theta + self.delta_theta_AB() + self.theta_misalign()

    def theta_C(self):  # todo
        pass

    def theta_open(self):
        return 2 * self.theta_misalign()


class Lengths:
    def __init__(self, chi, theta, nu, t, r, p):
        self.chi = chi
        self.theta = theta
        self.nu = nu
        self.t = t
        self.r = r
        self.p = p
        self.__angles = Angles(chi, theta, nu, t, r, p)

    def geo_focus(self):
        return np.cos(self.chi - self.theta) / (np.cos(self.chi + self.theta) / self.p + 2.0 / self.r)

    def single_ray_focus(self):
        return (self.r * np.sin(2.0 * self.theta)) / (
                2.0 * np.sin(self.chi + self.theta) + (1 + self.nu) * np.sin(2.0 * self.chi) * np.cos(
            self.chi + self.theta))

    def width(self):
        theta_open = self.__angles.theta_open()
        chi_B = self.__angles.chi_B()
        theta_B = self.__angles.theta_B()
        fg_B = Lengths(chi_B, theta_B, self.nu, self.t, self.r, self.p).geo_focus()
        return - fg_B * theta_open

    def foot_length(self):  # todo: double check using point B or point C for theta and chi
        width = self.width()
        theta_B = self.__angles.theta_B()
        chi_B = self.__angles.chi_B()
        return width / np.cos(theta_B - chi_B)


if __name__ == '__main__':
    # print(Angles(chi=np.radians(4.4671),theta=np.radians(8.99),nu=0.2,t=0.3,r=2000,p=22000).theta_misalign())
    import os

    print(os.path.dirname(os.path.abspath(__file__)))
