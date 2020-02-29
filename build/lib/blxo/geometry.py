import os

file_dir = os.path.dirname(os.path.abspath(__file__))
import sys

sys.path.extend([file_dir])
import numpy as np
from scipy.optimize import fsolve


class Angles:
    def __init__(self, chi, theta, nu, t, r, p):
        self.__chi = chi
        self.__theta = theta
        self.__nu = nu
        self.__t = t
        self.__r = r
        self.__p = p

    # magic condition core angles
    def delta_phi_AB(self):
        return -self.__t * np.tan(self.__theta - self.__chi) / self.__r

    def delta_chi_AB(self):
        return -(1 + self.__nu) * self.__t * np.sin(2 * self.__chi) / (2 * self.__r)

    def psi_AB(self):
        return (self.__t / self.__p) * (np.sin(2 * self.__theta) / np.cos(self.__chi - self.__theta))

    def delta_theta_AB(self):
        '''
        Notes:
            This $\Delta\theta_{AB}$ is caused by the finite source distance.
            ONLY when the 'magic condition' is met, this is the difference of Bragg angles between point A and B.
            When the 'magic condition' is not met, a misalignment angle should be added for the difference of Bragg
            angles between A and B.
        Returns:

        '''
        return -1 / 2 * self.psi_AB()

    # Other psi angles (angles from the finite source distance)
    def psi_BC(self):  # todo: double check using point A or point C for theta and chi
        qmb_foot_length = Lengths(self.__chi, self.__theta, self.__nu, self.__t, self.__r, self.__p).foot_length()
        return np.cos(self.theta_B() + self.chi_B()) * qmb_foot_length / self.__p

    def chi_B(self):
        return self.__chi + self.delta_chi_AB()

    def chi_C(self):
        return self.chi_B()

    def theta_misalign(self):
        try:
            from blxo.mc import magic_condition_angles
        except:
            from mc import magic_condition_angles
        return magic_condition_angles(self.__chi, self.__theta, self.__nu, self.__t, self.__r, self.__p)

    def theta_B(self):
        return self.__theta + self.delta_theta_AB() + self.theta_misalign()

    def theta_C(self):  # todo
        pass

    def theta_open(self):
        return 2 * self.theta_misalign()


class Lengths:
    def __init__(self, chi, theta, nu, t, r, p):
        self.__chi = chi
        self.__theta = theta
        self.__nu = nu
        self.__t = t
        self.__r = r
        self.__p = p
        self.__angles = Angles(chi, theta, nu, t, r, p)

    def geo_focus(self):
        return np.cos(self.__chi - self.__theta) / (np.cos(self.__chi + self.__theta) / self.__p + 2.0 / self.__r)

    def single_ray_focus(self):
        return (self.__r * np.sin(2.0 * self.__theta)) / (
                2.0 * np.sin(self.__chi + self.__theta) + (1 + self.__nu) * np.sin(2.0 * self.__chi) * np.cos(
            self.__chi + self.__theta))

    def width(self):
        theta_open = self.__angles.theta_open()
        chi_B = self.__angles.chi_B()
        theta_B = self.__angles.theta_B()
        fg_B = Lengths(chi_B, theta_B, self.__nu, self.__t, self.__r, self.__p).geo_focus()
        return - fg_B * theta_open

    def foot_length(self):  # todo: double check using point B or point C for theta and chi
        width = self.width()
        theta_B = self.__angles.theta_B()
        chi_B = self.__angles.chi_B()
        return width / np.cos(theta_B - chi_B)


if __name__ == '__main__':
    print(Angles(chi=np.radians(4.4671), theta=np.radians(8.99), nu=0.2, t=0.3, r=2000, p=22000).theta_misalign())
