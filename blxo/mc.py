import os

file_dir = os.path.dirname(os.path.abspath(__file__))
import sys

sys.path.extend([file_dir])
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

try:
    import blxo.geometry as geometry
except:
    import geometry

'''
psi: Incident ray angle divergence from a point source with FINITE source distance.
beta: Emergent ray opening angle caused by the misalignment from the magic condition. To be specific, 
        it is 2 times the Magic Condition Misalignment.
theta: Angles refering to Bragg angle.
delta_theta_AB: The Bragg angle divergence from A to B. Or `theta_B-theta_A`.
chi: Angles refering to asymmetry angle. 


TODO: - Bandwidth
      - Equation solvers.
      - Darwin width function should be a calculator instead of a dictionary reader.
'''


class BentLaueMono(object):
    def __init__(self, chi, theta, nu, t, r, p, s=0, hkl=(1, 1, 1)):  # radians, mm.
        self.__chi = chi
        self.__theta = theta
        self.__nu = nu
        self.__t = t
        self.__r = r
        self.__p = p
        self.__s = s
        self.__hkl = hkl
        self.__angles = geometry.Angles(chi, theta, nu, t, r, p)
        self.__angle_resolution = self.angle_resolution_function()
        self.__energy_resolution = self.energy_resolution_function()
        self.__lengths = geometry.Lengths(chi, theta, nu, t, r, p)
        self.__qmb = self.quasi_mono_beam()
        self.__focal_size = self.focal_size_function()
        self.__f2d_optimal = self.f2d_optimal_function()

    @property
    def qmb(self):
        return self.__qmb

    @qmb.setter
    def qmb(self, var):
        raise Exception('"qmb" value should not be changed directly.')

    @property
    def angles(self):
        return self.__angles

    @angles.setter
    def angles(self, var):
        raise Exception('"angles" value should not be changed directly.')

    @property
    def lengths(self):
        return self.__lengths

    @lengths.setter
    def lengths(self, var):
        raise Exception('"lengths" value should not be changed directly.')

    @property
    def angle_resolution(self):
        return self.__angle_resolution

    @angle_resolution.setter
    def angle_resolution(self, var):
        raise Exception('"angle_resolution" value should not be changed directly.')

    @property
    def energy_resolution(self):
        return self.__energy_resolution

    @energy_resolution.setter
    def energy_resolution(self, var):
        raise Exception('"energy_resolution" value should not be changed directly.')

    @property
    def focal_size(self):
        return self.__focal_size

    @focal_size.setter
    def focal_size(self, var):
        raise Exception('"focal_size" value should not be changed directly.')

    @property
    def f2d_optimal(self):
        return self.__f2d_optimal

    @f2d_optimal.setter
    def f2d_optimal(self, var):
        raise Exception('"f2d_optimal" value should not be changed directly.')

    def quasi_mono_beam(self):
        # lengths = geometry.Lengths(self.chi, self.theta, self.nu, self.t, self.r, self.p)
        width = self.__lengths.width()
        foot_length = self.__lengths.foot_length()

        def angular_spread(self):
            delta_theta_Q1 = self.angle_resolution['dtheta_1']
            delta_theta_Q2 = -(self.__angles.psi_AB() - self.__angles.psi_BC()) / 2
            delta_theta_Q = delta_theta_Q1 + delta_theta_Q2
            return delta_theta_Q

        ang_sprd = angular_spread(self)

        return {'width': width, 'foot_length': foot_length, 'angular_spread': ang_sprd}

    def angle_resolution_function(self):
        dtheta1 = -self.__t / self.__r * (np.cos(self.__chi) ** 2 - self.__nu * np.sin(self.__chi) ** 2) * np.tan(
            self.__theta) \
                  + (self.__chi - self.__chi) + (self.__theta - self.__theta) + (self.__nu - self.__nu) + (
                          self.__t - self.__t) + (self.__r - self.__r) + (self.__p - self.__p) + (
                          self.__s - self.__s)  # d-spacing
        dtheta2 = self.__angles.delta_theta_AB() \
                  + (self.__chi - self.__chi) + (self.__theta - self.__theta) + (self.__nu - self.__nu) + (
                          self.__t - self.__t) + (self.__r - self.__r) + (self.__p - self.__p) + (
                          self.__s - self.__s)  # finite source distance
        dtheta3 = darwin_width(self.__hkl) * np.tan(self.__theta) \
                  + (self.__chi - self.__chi) + (self.__r - self.__r) + (self.__nu - self.__nu) + (
                          self.__t - self.__t) + (self.__p - self.__p) + (self.__s - self.__s)  # darwin width
        dtheta4 = self.__s / self.__p \
                  + (self.__chi - self.__chi) + (self.__theta - self.__theta) + (self.__nu - self.__nu) + (
                          self.__t - self.__t) + (self.__r - self.__r)  # source size. source_size FWHM at BMIT_BM is 118um
        dtheta_all = np.sqrt(
            (dtheta1 + dtheta2) ** 2 + dtheta3 ** 2 + dtheta4 ** 2)  # Absolute value only. Non-directional.
        return {'dtheta_1': dtheta1,  # dtheta
                'dtheta_2': dtheta2,  # dtheta
                'dtheta_3': dtheta3,  # dtheta
                'dtheta_4': dtheta4,  # dtheta
                'dtheta_all': dtheta_all}

    def energy_resolution_function(self):
        angle_res = self.angle_resolution
        return {'de_1': abs(- angle_res['dtheta_1'] / np.tan(self.__theta)),  # dE/E
                'de_2': abs(- angle_res['dtheta_2'] / np.tan(self.__theta)),  # dE/E
                'de_3': abs(- angle_res['dtheta_3'] / np.tan(self.__theta)),  # dE/E
                'de_4': abs(- angle_res['dtheta_4'] / np.tan(self.__theta)),  # dE/E
                'de_all': angle_res['dtheta_all'] / np.tan(self.__theta)}  # Absolute value only. Non-directional.

    def f2d_optimal_function(self):
        # qmb = self.quasi_mono_beam()
        width = self.qmb['width']
        ang_sprd = self.qmb['angular_spread']
        fg = self.__lengths.geo_focus()
        f2d = -(1 / self.__r + np.cos(self.__chi + self.__theta) / self.__p) * width * fg / (
                ang_sprd * np.cos(self.__chi - self.__theta))
        return f2d

    def focal_size_function(self):
        fg = self.__lengths.geo_focus()
        width = self.qmb['width']
        return np.sqrt((fg * self.__s / self.__p) ** 2 + width ** 2)


def magic_condition_angles(chi, theta, nu, t, r, p):
    angles = geometry.Angles(chi, theta, nu, t, r, p)
    theta_misalignment = angles.delta_phi_AB() - angles.delta_chi_AB() + angles.delta_theta_AB()
    return theta_misalignment


def magic_condition_foci(chi, theta, nu, t, r, p):
    lengths = geometry.Lengths(chi, theta, nu, t, r, p)
    return lengths.geo_focus() - lengths.single_ray_focus()


def darwin_width(hkl=(1, 1, 1)):
    darwin_code = {(1, 1, 1): 134.3e-6,
                   (2,2,0): 56.62,
                   (3,1,1): 26.71,
                   (4,0,0): 23.72,
                   (3,3,1): 13.59,
                   (4,2,2): 13.85,
                   (3,3,3): 8.25,
                   (5,1,1): 8.25,
                   (4,4,0): 9.10,
                   (5,3,1): 5.60,
                   (6,2,0): 6.41,
                   (5,3,3): 4.02,
                   (4,4,4): 4.74,
                   (5,5,1): 3.01,
                   (7,1,1): 3.01,
                   (6,4,2): 3.62,
                   (5,5,3): 2.33,
                   (7,3,1): 2.33,
                   (8,0,0): 2.83,
                   (7,3,3): 1.84,
                   (6,6,0): 2.27,
                   (8,2,2): 2.27,
                   (5,5,5): 1.48,
                   (7,5,1): 1.48,
                   (8,4,0): 1.84,
                   (7,5,3): 1.22,
                   (9,1,1): 1.22,
                   (6,6,4): 1.53,
                   (9,3,1): 1.03,
                   (8,4,4): 1.28,
                   (7,5,5): 0.86,
                   (9,3,3): 0.86,
                   (7,7,1): 0.86,
                   (8,6,2): 1.09,
                   (10,2,0): 1.09,
                   (9,5,1): 0.72,
                   (7,7,3): 0.72,
                   (9,5,3): 0.61,
                   (10,4,2): 0.78,
                   (11,1,1): 0.53,
                   (7,7,5): 0.53,
                   (8,8,0): 0.70}
    hkl_sum = sum([i**2 for i in hkl])
    for key in darwin_code:
        if sum([i**2 for i in key]) == hkl_sum:
            return darwin_code[key]
    raise ValueError('Reflection ' + str(hkl) + ' is not valid.')
    return


if __name__ == '__main__':
    print(geometry.Angles(chi=np.radians(10), theta=np.radians(8.99), nu=0.2, t=0.3, r=2000,
                          p=22000).theta_misalign())
    mono = BentLaueMono(chi=np.radians(10), theta=np.radians(8.99), nu=0.2, t=0.3, r=2000, p=22000,hkl=(1,2,2))
    print(mono.quasi_mono_beam()['width'])
    print(mono.quasi_mono_beam()['angular_spread'])
