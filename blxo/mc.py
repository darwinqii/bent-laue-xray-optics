# %%
import os

file_dir = os.path.dirname(os.path.abspath(__file__))
import sys

sys.path.extend([file_dir])
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
import blxo.geometry as geometry

'''
alpha: Incident ray angle divergence from a point source with FINITE source distance.
beta: Emergent ray opening angle caused by the misalignment from the magic condition. To be specific, 
        it is 2 times the Magic Condition Misalignment.
theta: Angles refering to Bragg angle.
delta_theta_AB: The Bragg angle divergence from A to B. Or `theta_B-theta_A`.
chi: Angles refering to asymmetry angle. 


TODO: keep the calculated values (eg. geo_focus()) as a method, so that it cannot be altered simply by assigning a value

'''


class BentLaueMono(object):
    chi = 5

    def __init__(self, chi, theta, nu, t, r, p, s=0, hkl=(1, 1, 1)):  # radians, mm.
        self.chi = chi
        self.theta = theta
        self.nu = nu
        self.t = t
        self.r = r
        self.p = p
        self.s = s
        self.hkl = hkl
        self.__angles = geometry.Angles(chi, theta, nu, t, r, p)
        self.__lengths = geometry.Lengths(chi, theta, nu, t, r, p)
        self.qmb = self.quasi_mono_beam()
        # self.angle_resolution=None
        # self.energy_resolution=None

        # self.geo_focus = geo_focus(chi=chi, theta=theta, p=p, r=r)
        # self.single_ray_focus = single_ray_focus(chi, theta, nu, r)
        # self.theta_Bmis = magic_condition_angles(chi, theta, nu, t, r, p)
        # self.energy_resolution = self.EnergyResolution(self)
        # self.quasi_mono_beam = QuasiMonoBeam(chi, theta, nu, t, r, p)

    def quasi_mono_beam(self):
        # lengths = geometry.Lengths(self.chi, self.theta, self.nu, self.t, self.r, self.p)
        width = self.__lengths.width()
        foot_length = self.__lengths.foot_length()

        def angular_spread(self):
            delta_theta_Q1 = self.angle_resolution()['dtheta_1']
            delta_theta_Q2 = -(self.__angles.alpha_AB() - self.__angles.alpha_BC()) / 2
            delta_theta_Q = delta_theta_Q1 + delta_theta_Q2
            return delta_theta_Q

        ang_sprd = angular_spread(self)

        return {'width': width, 'foot_length': foot_length, 'angular_spread': ang_sprd}

    def angle_resolution(self):
        dtheta1 = -self.t / self.r * (np.cos(self.chi) ** 2 - self.nu * np.sin(self.chi) ** 2) * np.tan(
            self.theta)  # d-spacing
        dtheta2 = self.__angles.delta_theta_AB()  # finite source distance
        dtheta3 = darwin_width(self.hkl) * np.tan(self.theta)  # darwin width
        dtheta4 = self.s / self.p  # source size
        dtheta_all = np.sqrt(
            (dtheta1 + dtheta2) ** 2 + dtheta3 ** 2 + dtheta4 ** 2)  # Absolute value only. Non-directional.
        return {'dtheta_1': dtheta1,
                'dtheta_2': dtheta2,
                'dtheta_3': dtheta3,
                'dtheta_4': dtheta4,
                'dtheta_all': dtheta_all}

    def energy_resolution(self):
        angle_res = self.angle_resoltution()
        return {'de_1': - angle_res['dtheta_1'] / np.tan(self.theta),
                'de_2': - angle_res['dtheta_2'] / np.tan(self.theta),
                'de_3': - angle_res['dtheta_3'] / np.tan(self.theta),
                'de_4': - angle_res['dtheta_4'] / np.tan(self.theta),
                'de_all': angle_res['dtheta_all'] / np.tan(self.theta)}  # Absolute value only. Non-directional.

    def f2d_optimal(self):
        # qmb = self.quasi_mono_beam()
        width = self.qmb['width']
        ang_sprd = self.qmb['angular_spread']
        fg = self.__lengths.geo_focus()
        f2d = -(1 / self.r + np.cos(self.chi + self.theta) / self.p) * width * fg / (
                ang_sprd * np.cos(self.chi - self.theta))
        return f2d

    def focal_size(self):
        fg = self.__lengths.geo_focus()
        width = self.qmb['width']
        return np.sqrt((fg * self.s / self.p) ** 2 + width ** 2)


def magic_condition_angles(chi, theta, nu, t, r, p):
    print('cool')
    angles = geometry.Angles(chi, theta, nu, t, r, p)
    theta_misalignment = angles.delta_psi_AB() - angles.delta_chi_AB() + angles.delta_theta_AB()
    return theta_misalignment


def darwin_width(hkl=(1, 1, 1)):
    darwin_code = {(1, 1, 1): 134.3e-6}
    return darwin_code[hkl]


if __name__ == '__main__':
    print(geometry.Angles(chi=np.radians(10), theta=np.radians(8.99), nu=0.2, t=0.3, r=2000,
                          p=22000).theta_misalign())
    mono = BentLaueMono(chi=np.radians(10), theta=np.radians(8.99), nu=0.2, t=0.3, r=2000, p=22000)
    print(mono.quasi_mono_beam()['width'])
    print(mono.quasi_mono_beam()['angular_spread'])


    #
    # def magic_condition_foci(chi, theta, nu, r, p):
    #     return single_ray_focus(chi, theta, nu, r) - geo_focus(chi, theta, p, r)

    #
    #
    # def theta_B2_fsolver(chi, theta, nu, t, r, p, verbose=False):  # TODO: double check this function
    #     def theta_B2_equation(theta_b2, chi=chi, theta=theta, nu=nu, r=r, p=p, t=t, verbose=verbose):
    #         # print(chi)
    #         theta_offset = magic_condition_angles(chi, theta, nu, t, r, p)
    #         theta_sd = -(t / (2 * p)) * (np.sin(2 * theta) / np.cos(chi - theta))
    #         theta_b1 = theta + theta_sd + theta_offset
    #         delta_chi = -(t / r) * ((1 + nu) / 2) * np.sin(2 * chi)
    #         chi1 = chi + delta_chi
    #         fg1 = geo_focus(chi1, theta_b1, r, p)
    #         length = -fg1 * 2 * theta_offset / np.cos(theta_b2 - chi1)
    #         equation = theta_b2 - theta_b1 - length * np.cos(chi1 + theta_b2) / (2 * p)
    #         # print(t)
    #         if verbose:
    #             print('\noffset:', theta_offset)
    #             print('theta_sd:', theta_sd)
    #             print('theta_b0:', theta)
    #             print('theta_b1:', theta_b1)
    #             print('delta_chi:', delta_chi)
    #             print('chi1:', chi1)
    #             print('arc length:', length)
    #             print('theta_b2', theta_b2)
    #         return equation
    #
    #     for arg in [chi, theta, nu, r, p, t]:
    #         if not np.isscalar(arg):  # if arg is a sequence (list or array), theta will be expanded to the same
    #         length
    #             arg_length = len(arg)
    #             break
    #         else:
    #             arg_length = 1
    #     theta_B2 = fsolve(theta_B2_equation, np.ones(arg_length) * theta, args=(chi, theta, nu, r, p, t))[0]
    #     return theta_B2
    #
    #
    # def focal_size(chi, theta, nu, T, R, D, S0=0, verbose=False):
    #     f_g = geo_focus(chi, theta, R, D)
    #     width = mono_beam_width(chi, theta, nu, T, R, D)
    #     size = abs((f_g / D) * S0) + abs(width)
    #     return size
    #
    #
    # def energy_resolution_all(chi, theta, nu, R, D, T,
    #                           S=0.118,  # source_size FWHM at BMIT_BM is 118um
    #                           hkl=[1, 1, 1],  # Darwin width for silicon [111] is dE/E=134.3 x 10^-6
    #                           return_all=False):
    #     # energy resolution contributions
    #     term1 = delta_theta_d_spacing(chi, theta, nu, R, T) / (-np.tan(theta))  # d-spacing
    #     term2 = delta_theta_source_distance(chi, theta, D, T) / (-np.tan(theta))  # finite source distance
    #     term3 = delta_theta_darwin(theta, hkl) / (-np.tan(theta))  # darwin width
    #     term4 = delta_theta_source_size(D, S) / (-np.tan(theta))  # source size
    #
    #     resolution = np.sqrt((term1 + term2) ** 2 + term3 ** 2 + term4 ** 2)
    #
    #     if return_all:
    #         return resolution, term1, term2, term3, term4
    #     else:
    #         return resolution
    #
    #
    # def mc_angle_solver(theta=None, chi=None, R=None, nu=None, D=None, T=None, verbose=False):
    #     chi_range_a = np.radians(np.arange(-45.0, 45.0, 0.001))
    #     scores_range_a = magic_condition_angles(chi_range_a, theta=theta, T=T, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    #     # print(np.degrees(ulti_chi))
    #     if verbose:
    #         plt.plot(np.degrees(chi_range_a), scores_range_a)
    #         # plt.ylim(-0.0005,0.0005)
    #         plt.grid()
    #
    #     chi_range_b = np.arange(ulti_chi - np.radians(1), ulti_chi + np.radians(1), 0.00001)
    #     scores_range_b = magic_condition_angles(chi_range_b, theta=theta, T=T, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    #
    #     chi_range_c = np.arange(ulti_chi - np.radians(0.0001), ulti_chi + np.radians(0.0001), 0.0000001)
    #     scores_range_c = magic_condition_angles(chi_range_c, theta=theta, T=T, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    #
    #     return ulti_chi
    #
    #
    # def mc_focus_solver(theta=None, chi=None, R=None, nu=None, D=None, T=None):
    #     chi_range_a = np.radians(np.arange(-80.0, 80.0, 0.01))
    #     scores_range_a = magic_condition_foci(chi_range_a, theta=theta, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    #     plt.figure()
    #     plt.plot(chi_range_a, scores_range_a)
    #     plt.ylim(-10000000, 10000000)
    #     plt.grid()
    #
    #     chi_range_b = np.arange(ulti_chi - np.radians(1), ulti_chi + np.radians(1), 0.0001)
    #     scores_range_b = magic_condition_foci(chi_range_b, theta=theta, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    #
    #     chi_range_c = np.arange(ulti_chi - np.radians(0.0001), ulti_chi + np.radians(0.0001), 0.0000001)
    #     scores_range_c = magic_condition_foci(chi_range_c, theta=theta, D=D, R=R, nu=nu)
    #     ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    #     return ulti_chi
    #
    #
    # def radius_from_geo_focus(f_g, chi, theta, D):
    #     chi = np.radians(chi)
    #     theta = np.radians(theta)
    #     R = 2 / (np.cos(chi - theta) / f_g - np.cos(chi + theta) / D)
    #     return R
    #
    #
    # def mc_solver_chi_focus(chi=None, f_g=None, theta=None, nu=None, D=None, T=None):
    #     ''' Given f_g (instead of R), solve the chi angle for magic condition
    #     '''
    #     chi_range_a = np.arange(-90.0, 90.0, 0.01)
    #     R = radius_from_geo_focus(f_g, chi_range_a, theta=theta, D=D)
    #     scores_range_a = magic_condition_angles(chi_range_a, theta=theta, T=T, D=D, R=R, nu=nu)
    #     magic_chi = chi_range_a[abs(scores_range_a).argmin()]
    #     #     plt.plot(chi_range_a,scores_range_a)
    #     #     plt.ylim(-0.0005,0.0005)
    #     #     plt.grid()
    #
    #     chi_range_b = np.arange(magic_chi - 1, magic_chi + 1, 0.0001)
    #     R = radius_from_geo_focus(f_g, chi_range_b, theta=theta, D=D)
    #     scores_range_b = magic_condition_angles(chi_range_b, theta=theta, T=T, D=D, R=R, nu=nu)
    #     magic_chi = chi_range_b[abs(scores_range_b).argmin()]
    #
    #     chi_range_c = np.arange(magic_chi - 0.0001, magic_chi + 0.0001, 0.0000001)
    #     R = radius_from_geo_focus(f_g, chi_range_c, theta=theta, D=D)
    #     scores_range_c = magic_condition_angles(chi_range_c, theta=theta, T=T, D=D, R=R, nu=nu)
    #     magic_chi = chi_range_c[abs(scores_range_c).argmin()]
    #
    #     return magic_chi
    #
    #
    # def magic_condition_bend_radius(chi, theta, D, nu):
    #     R = ((1 + nu) * D * np.sin(2 * chi) * np.cos(chi - theta) - 2 * D * np.sin(theta - chi)) / np.sin(2 * theta)
    #     return R
    #
    #
    # def magic_condition_source_distance(chi, theta, R, nu):
    #     D = R * np.sin(2 * theta) / ((1 + nu) * np.sin(2 * chi) * np.cos(chi - theta) + 2 * np.sin(chi - theta))
    #     return D
    #
    #
    # def energy_resolution(chi, theta, R, D, T, nu):
    #     theta = np.radians(theta)
    #     chi = np.radians(chi)
    #     term1 = T * (np.cos(chi) ** 2 - nu * np.sin(chi) ** 2) / R
    #     term2 = T * (np.cos(theta) ** 2) / (2 * D * np.cos(chi - theta))
    #     resolution = term1 + term2
    #     return resolution
    #
    #
    # def ulti_magi_condi(chi, theta=8.99, nu=0.2, T=0.3, D=2200, return_all=False):
    #     theta = np.radians(theta)
    #     chi = np.radians(chi)
    #     R = ((1 + nu) * D * np.sin(2 * chi) * np.cos(chi - theta) - 2 * D * np.sin(theta - chi)) / np.sin(2 * theta)
    #     try:
    #         len(R)
    #     except:
    #         print('R:', R)
    #     # energy resolution equation
    #     term1 = T * (np.cos(chi) ** 2 - nu * np.sin(chi) ** 2) / R
    #     term2 = T * (np.cos(theta) ** 2) / (2 * D * np.cos(chi - theta))
    #     resolution = term1 + term2
    #     if return_all:
    #         return resolution, term1, term2
    #     else:  # so that we can use this function with `fsolve`
    #         return resolution
    #
    #
    # def ulti_magi_condi_with_R(chi, theta=8.99, nu=0.2, T=0.3, R=2000, return_all=False):
    #     theta = np.radians(theta)
    #     chi = np.radians(chi)
    #     D = R * np.sin(2 * theta) / ((1 + nu) * np.sin(2 * chi) * np.cos(chi - theta) + 2 * np.sin(chi - theta))
    #     try:
    #         len(D)
    #     except:
    #         print('D:', D)
    #     # energy resolution equation
    #     term1 = T * (np.cos(chi) ** 2 - nu * np.sin(chi) ** 2) / R
    #     term2 = T * (np.cos(theta) ** 2) / (2 * D * np.cos(chi - theta))
    #     #     term3 = term1*term2
    #     resolution = term1 + term2
    #     if return_all:
    #         return resolution, term1, term2
    #     else:  # so that we can use this function with `fsolve`
    #         return resolution
    #
    #
    # def ulti_magi_condi_chi_solver(theta=8.98, T=0.3, D=22000, nu=0.2):
    #     chi_range_a = np.arange(-45.0, 45.0, 0.01)
    #     scores_range_a = ulti_magi_condi(chi_range_a, theta=theta, T=T, D=D, nu=nu)
    #     ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    #     plt.plot(chi_range_a, scores_range_a)
    #     plt.ylim(-0.0001, 0.0001)
    #     plt.ylabel('$\Delta E/E$')
    #     plt.grid()
    #     chi_range_b = np.arange(ulti_chi - 1, ulti_chi + 1, 0.0001)
    #     scores_range_b = ulti_magi_condi(chi_range_b, theta=theta, T=T, D=D, nu=nu)
    #     ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    #
    #     #     chi_range_c = np.arange(ulti_chi-0.0001,ulti_chi+0.0001,0.0000001)
    #     #     scores_range_c=ulti_magi_condi(chi_range_c,theta=theta,T=T,D=D,nu=nu)
    #     #     ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    #
    #     return ulti_chi
    #
    #
    # def energy_divergence(theta1, theta2, hkl=[1, 1, 1]):
    #     import math_physics as mp
    #     energy1 = mp.bragg(hkl=hkl, theta=theta1)[0]
    #     energy2 = mp.bragg(hkl=hkl, theta=theta2)[0]
    #     energy_divergence = (energy2 - energy1) * 1000  # return in eV
    #     return energy_divergence
    #
    #
    # def dy_to_theta(dy, theta, chi, R, D=np.Infinity):
    #     theta = np.radians(theta)
    #     chi = np.radians(chi)
    #     y0 = -np.sin(theta + chi) * R
    #     theta_dy = np.arcsin((-y0 - dy) / R) - chi - np.arcsin(dy / D)
    #     theta_dy = np.degrees(theta_dy)
    #     return theta_dy
    #
    #
    # def mc_y_evaluation(chi=4.4671, theta=8.99, nu=0.2, T=0.3, R=2000, D=22000, d_y=np.arange(-2.5, 2.5, 0.001),
    #                     verbose=False):
    #     theta_0 = np.radians(theta)
    #     chi = np.radians(chi)
    #     y_0 = -np.sin(theta_0 + chi) * R
    #     print(y_0)
    #
    #     d_y = d_y
    #     theta = np.arcsin((-y_0 - d_y) / R) - chi - np.arcsin(d_y / D)  # radians
    #     mc_angle_distance = magic_condition_angle(np.degrees(chi), np.degrees(theta), nu, T, R, D)  # in radians
    #     mbw = mono_beam_width(np.degrees(chi), np.degrees(theta), nu, T, R, D)
    #     if verbose:
    #         # Plot
    #         plt.figure(figsize=(13, 5))
    #         plt.subplot(1, 2, 1)
    #         plt.plot(d_y, np.degrees(mc_angle_distance), label='Opening angle at Exit')
    #         plt.legend()
    #         plt.ylabel('deg')
    #         plt.xlabel('Relative vertical position (mm)')
    #         plt.title('Opening angle vs Vertical x-ray position')
    #         plt.subplot(1, 2, 2)
    #         plt.plot(d_y, mbw, color='r', label='Monochromatic beam width (mm)')
    #         plt.xlabel('Relative vertical position (mm)')
    #         plt.legend()
    #         plt.title('Monochromatic beam width vs Vertical x-ray position')
    #
    #     # Theta angle divergence over the whole beam vertically
    #     try:
    #         D_theta = np.degrees(theta[-1] - theta[0])  # in degree
    #
    #         if verbose:
    #             print('Incident theta angle divergence over the whole beam vertically:\n', abs(D_theta))
    #             print('The biggest exiting point angle divergence is (deg):\n',
    #                   np.degrees(max(abs(mc_angle_distance[-1]), abs(mc_angle_distance[0]))))
    #             D_theta_2 = 5 / (R * np.cos(chi + theta_0)) + 5 / D
    #             print(np.degrees(D_theta_2))
    #         #             print(np.radians(D_theta)/np.tan(theta_0)*12.658)
    #         D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0] - mp.bragg(theta=np.degrees(theta[0]))[0]
    #         if verbose:
    #             print('Overall energy range:\n', D_energy)
    #     except:
    #         pass
    #
    #     return mc_angle_distance, mbw
    #
    #
    # def mc_theta_evaluation(chi=4.4671, theta=8.99, nu=0.2, T=0.3, R=2000, D=22000, verbose=False):
    #     theta_0 = np.radians(theta)
    #     chi = np.radians(chi)
    #
    #     d_theta = np.arange(-np.pi / 4, np.pi / 4, 0.0001)
    #     theta = theta_0 + d_theta  # in radians
    #     mc_angle_distance = magic_condition_angle(np.degrees(chi), np.degrees(theta), nu, T, R, D)  # in radians
    #     mbw = mono_beam_width(np.degrees(chi), np.degrees(theta), nu, T, R, D)
    #     if verbose:
    #         # Plot
    #         #         plt.figure(figsize=(13,5))
    #         #         plt.subplot(1,2,1)
    #         #         plt.plot(np.degrees(d_theta),np.degrees(mc_angle_distance),label='Opening angle at Exit')
    #         #         plt.legend()
    #         #         plt.ylabel('deg')
    #         #         plt.xlabel('Angle distance from Magic Condition Bragg angle (deg)')
    #         #         plt.title('Opening angle vs Bragg angle')
    #         #         plt.subplot(1,2,2)
    #         plt.plot(np.degrees(d_theta), mbw * 1000, color='r', label='Monochromatic beam width ($\mu m$)')
    #         plt.xlabel(r'Relative distance from magic condition Bragg angle $\theta_B$ (deg)')
    #         plt.ylabel('Monochromatic beam width ($\mu m$)')
    #         #         plt.legend()
    #         plt.title('Monochromatic beam width vs Bragg angle')
    #         plt.grid()
    #
    #     #     # Theta angle divergence over the whole beam vertically
    #     #     D_theta = np.degrees(d_theta[-1]-d_theta[0]) # in degree
    #     #     if verbose:
    #     #         print('Incident theta angle divergence over the whole beam vertically:\n',abs(D_theta))
    #     #         print('The biggest exiting point angle divergence is (deg):\n', np.degrees(max(abs(
    #     mc_angle_distance[
    #     #         -1]),abs(mc_angle_distance[0]))))
    #
    #     #     D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0]-mp.bragg(theta=np.degrees(theta[0]))[0]
    #     # #     print(theta[-1],theta[0])
    #     #     if verbose:
    #     #         print('Overall energy range: (keV)\n',D_energy)
    #
    #     return mc_angle_distance, None, None
    #
    #
    # def mc_chi_evaluation(chi=4.4671, theta=8.99, nu=0.2, T=0.3, R=2000, D=22000, verbose=False):
    #     theta = np.radians(theta)
    #     chi_0 = np.radians(chi)
    #
    #     d_chi = np.arange(-np.pi / 4, np.pi / 4, 0.0001)
    #     chi = chi_0 + d_chi  # in radians
    #     mc_angle_distance = magic_condition_angle(np.degrees(chi), np.degrees(theta), nu, T, R, D)  # in radians
    #     mbw = mono_beam_width(np.degrees(chi), np.degrees(theta), nu, T, R, D)
    #     if verbose:
    #         #         # Plot
    #         #         plt.figure(figsize=(13,5))
    #         #         plt.subplot(1,2,1)
    #         #         plt.plot(np.degrees(d_chi),np.degrees(mc_angle_distance),label='Opening angle at Exit')
    #         #         plt.legend()
    #         #         plt.ylabel('deg')
    #         #         plt.xlabel('Angle distance from Magic Condition Chi angle (deg)')
    #         #         plt.title('Opening angle VS asymmetry angle')
    #         #         plt.subplot(1,2,2)
    #         plt.plot(np.degrees(d_chi), mbw * 1000, color='r')
    #         plt.xlabel('Relative distance from magic condition asymmetry angle $\chi$ (deg)')
    #         plt.ylabel('Monochromatic beam width ($\mu m$)')
    #         #         plt.legend()
    #         plt.title('Monochromatic beam width VS asymmetry angle')
    #         plt.grid()
    #
    #     #     # Theta angle divergence over the whole beam vertically
    #     #     D_theta = np.degrees(d_theta[-1]-d_theta[0]) # in degree
    #     #     if verbose:
    #     #         print('Incident theta angle divergence over the whole beam vertically:\n',abs(D_theta))
    #     #         print('The biggest exiting point angle divergence is (deg):\n', np.degrees(max(abs(
    #     mc_angle_distance[
    #     #         -1]),abs(mc_angle_distance[0]))))
    #
    #     #     D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0]-mp.bragg(theta=np.degrees(theta[0]))[0]
    #     #     print(theta[-1],theta[0])
    #     #     if verbose:
    #     #         print('Overall energy range: (keV)\n',D_energy)
    #
    #     return mc_angle_distance, None, None
    #
    #
    # def mc_theta_wanderoff(chi=4.4671, theta=8.99, nu=0.2, T=0.3, R=2000, D=22000, verbose=False):
    #     theta_0 = np.radians(theta)
    #     chi = np.radians(chi)
    #     d_theta = np.arange(-np.pi / 4, np.pi / 4, 0.0001)
    #     theta = theta_0 + d_theta  # in radians
    #     mc_angle_distance = magic_condition_angles(np.degrees(chi), np.degrees(theta), nu, T, R, D)  # in radians
    #     mbw = mono_beam_width(np.degrees(chi), np.degrees(theta), nu, T, R, D)
    #     return np.degrees(d_theta), mbw * 1000
    #
    #
    # def mc_chi_wanderoff(chi=4.4671, theta=8.99, nu=0.2, T=0.3, R=2000, D=22000, verbose=False):
    #     theta = np.radians(theta)
    #     chi_0 = np.radians(chi)
    #     d_chi = np.arange(-np.pi / 4, np.pi / 4, 0.0001)
    #     chi = chi_0 + d_chi  # in radians
    #     mc_angle_distance = magic_condition_angles(np.degrees(chi), np.degrees(theta), nu, T, R, D)  # in radians
    #     mbw = mono_beam_width(np.degrees(chi), np.degrees(theta), nu, T, R, D)
    #     return np.degrees(d_chi), mbw * 1000
