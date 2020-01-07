import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


# For the Magic Condition equation (angles)
def delta_chi(chi,nu,R,T,verbose=False):
    d_chi = -(1+nu)*T*np.sin(2*chi)/(2*R)
    if verbose:
        print('$\Delta\chi$:',d_chi)
    return d_chi

def delta_phi(chi,theta,R,T,verbose=False):
    d_phi = T*np.tan(theta-chi)/R
    if verbose:
        print('$\Delta\phi$:',d_phi)
    return d_phi

def delta_theta_source_distance(chi,theta,D,T,verbose=False):
    d_theta=-(T/(2*D))*(np.sin(2*theta)/np.cos(chi-theta))
    if verbose:
        print('Angular spread from source distance (Bragg angle variation):',d_theta)
    return d_theta

def delta_theta_d_spacing(chi,theta,nu,R,T,verbose=False):
    d_theta = -(T/R)*(np.cos(chi)**2-nu*np.sin(chi)**2)*np.tan(theta)
    if verbose:
        print('Angular spread from d-spacing',d_theta)
    return d_theta

def delta_theta_darwin(theta,hkl=[1,1,1],verbose=False):
    darwin_code = {[1,1,1]:134.3e-6}
    d_theta = darwin_code[hkl]*np.tan(theta)
    if verbose:
        print('Angular spread from Darwin width',d_theta)
    return d_theta

def delta_theta_source_size(D,S,verbose=False):
    d_theta = S/D
    if verbose:
        print('Angular spread from source size',d_theta)
    return d_theta

def geo_focus(chi,theta,R,D,verbose=False):
    f_g = np.cos(chi-theta)/(np.cos(chi+theta)/D+2.0/R)
    if verbose:
        print('Geometric focus (f_g)',f_g,'(mm)')
    return f_g

def single_ray_focus(chi,theta,R,D,nu,verbose=False):
    f_s = (R*np.sin(2.0*theta))/(2.0*np.sin(chi+theta)+(1+nu)*np.sin(2.0*chi)*np.cos(chi+theta))
    if verbose:
        print('Single ray focus (f_s)',f_s,'(mm)')
    return f_s

def magic_condition_angles(chi,theta,nu,T,R,D,verbose=False): # chi: rad; theta: rad; T,R,D: mm
    # todo: Correction for magic condition sign
    term1 = delta_chi(chi,nu,R,T,verbose)
    term2 = delta_phi(chi,theta,R,T,verbose)
    term3 = -delta_theta_source_distance(chi,theta,D,T,verbose)
    angles_distance = (term1+term2+term3)
    # print('I am working')
    if verbose:
        print('Contributions from chi, phi, source distance:',term1,term2,term3)
        print('Magic condition offset in angle:',angles_distance)
    return angles_distance

def magic_condition_foci(chi, theta,nu, R, D,verbose=False):
    # polychromatic focus function
    f_s = single_ray_focus(chi,theta,R,D,nu,verbose)
    # geometric focus function
    f_g = geo_focus(chi,theta,R,D,verbose)
    # magic condition
    foci_distance = f_s-f_g
    if verbose:
        print('Magic condition offset in foci distance:',foci_distance)
    return foci_distance

def theta_b2_fsolver(chi,theta,nu,T,R,D,verbose=False):
    def theta_b2_equation(theta_b2,chi=chi,theta=theta,nu=nu,R=R,D=D,T=T,verbose=verbose):
        # print(chi)
        theta_offset = magic_condition_angles(chi,theta,nu,T,R,D,verbose)
        theta_sd = -(T/(2*D))*(np.sin(2*theta)/np.cos(chi-theta))
        theta_b1 = theta+theta_offset+theta_sd
        delta_chi = -(T/R)*((1+nu)/2)*np.sin(2*chi)
        chi1 = chi+delta_chi
        fg1=geo_focus(chi1,theta_b1,R,D,verbose)
        length = -fg1*2*theta_offset/np.cos(theta_b2-chi1)
        equation = theta_b2-theta_b1-length*np.cos(chi1+theta_b2)/(2*D)
        if verbose:
            print('\noffset:',theta_offset)
            print('theta_sd:',theta_sd)
            print('theta_b0:',theta)
            print('theta_b1:',theta_b1)
            print('delta_chi:',delta_chi)
            print('chi1:',chi1)
            print('arc length:',length)
            print('theta_b2',theta_b2)
        return equation
    for arg in [chi,theta,nu,R,D,T]:
        if not np.isscalar(arg): # if arg is a sequence (list or array), theta will be expanded to the same length
            arg_length=len(arg)
            break
        else: arg_length=1
    theta_b2 = fsolve(theta_b2_equation,np.ones(arg_length)*theta,args=(chi,theta,nu,R,D,T))
    return theta_b2

def mono_beam_width(chi,theta,nu,T,R,D,verbose=False):
    # width = mono_beam_arc_length(chi,theta,nu,T,R,D)*np.cos(theta-chi)
    # todo: Correction for magic condition sign
    d_chi = delta_chi(chi,nu,R,T,verbose)
    chi1=chi+d_chi
    theta_offset = magic_condition_angles(chi,theta,nu,T,R,D,verbose)
    theta_sd = delta_theta_source_distance(chi,theta,D,T,verbose)
    theta_b1 = theta+theta_offset+theta_sd
    fg1=geo_focus(chi1,theta_b1,R,D,verbose)
    theta_open= 2*magic_condition_angles(chi,theta,nu,T,R,D,verbose)
    width = -fg1*theta_open
    if verbose:
        print('Quasi-mono beam width:',width)
    return width

def mono_beam_arc_length(chi,theta,nu,T,R,D,verbose=False):
    '''
    The length of the footprint of the quasi-mono beam on the exiting surface of the crystal.
    '''
    width = mono_beam_width(chi,theta,nu,T,R,D,verbose)
    chi2 = chi+delta_chi(chi,nu,R,T,verbose)
    theta_b2 = theta_b2_fsolver(chi,theta,nu,T,R,D,verbose)
    length = width/np.cos(theta_b2-chi2)
    if verbose:
        print('Quasi-mono beam footprint length:',length)
    return length

def mono_beam_energy_spread(chi,theta,nu,T,R,D,verbose=False):
    def mono_beam_delta_theta_source_distance(chi,theta,nu,T,R,D,verbose=False):
        length = mono_beam_arc_length(chi,theta,nu,T,R,D,verbose)
        delta_theta_B = -1/2 *(-2*delta_theta_source_distance(chi,theta,D,T,verbose)-length*np.cos(chi+theta)/D)
        if verbose:
            print('Quasi-mono beam Bragg angle variation caused by finite source distance',delta_theta_B)
        return delta_theta_B

    theta_d_spacing = delta_theta_d_spacing(chi,theta,nu,R,T,verbose)
    theta_source_distance=mono_beam_delta_theta_source_distance(chi,theta,nu,T,R,D,verbose)
    energy_spread = theta_d_spacing+theta_source_distance
    if verbose:
        print('Quasi-mono beam angular energy spread:',energy_spread)
    return energy_spread

def detector_distance_optimal(chi,theta,nu,T,R,D,verbose=False):
    width = mono_beam_width(chi,theta,nu,T,R,D,verbose)
    energy_spread = mono_beam_energy_spread(chi,theta,nu,T,R,D,verbose)
    fg=geo_focus(chi,theta,R,D,verbose)

    detector_distance = -(1/R+np.cos(chi+theta)/D)*width*fg/(energy_spread*np.cos(chi-theta))
    if verbose:
        print('Optimal detector distance:',detector_distance)
    return detector_distance

def focal_size(chi,theta,nu,T,R,D,S0=0,verbose=False):
    f_g = geo_focus(chi,theta,R,D)
    width = mono_beam_width(chi,theta,nu,T,R,D)
    size = abs((f_g/D)*S0) + abs(width)
    return size

def energy_resolution_all(chi,theta,nu,R,D,T,
                            S=0.118, # source_size FWHM at BMIT_BM is 118um
                            hkl=[1,1,1], # Darwin width for silicon [111] is dE/E=134.3 x 10^-6
                            return_all=False):
    # energy resolution contributions
    term1 = delta_theta_d_spacing(chi,theta,nu,R,T)/(-np.tan(theta))    # d-spacing
    term2 = delta_theta_source_distance(chi,theta,D,T)/(-np.tan(theta)) # finite source distance
    term3 = delta_theta_darwin(theta,hkl)/(-np.tan(theta))              # darwin width
    term4 = delta_theta_source_size(D,S)/(-np.tan(theta))               # source size
    
    resolution = np.sqrt((term1+term2)**2+term3**2+term4**2)
    
    if return_all:
        return resolution, term1, term2, term3, term4
    else:
        return resolution

def mc_angle_solver(theta=None, chi = None, R=None, nu=None,D=None,T=None,verbose=False):
    
    chi_range_a = np.radians(np.arange(-45.0,45.0,0.001))
    scores_range_a =  magic_condition_angles(chi_range_a,theta=theta,T=T,D=D,R=R,nu=nu)
    ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    # print(np.degrees(ulti_chi))
    if verbose:
        plt.plot(np.degrees(chi_range_a),scores_range_a)
        # plt.ylim(-0.0005,0.0005)
        plt.grid()
    
    chi_range_b = np.arange(ulti_chi-np.radians(1),ulti_chi+np.radians(1),0.00001)
    scores_range_b=magic_condition_angles(chi_range_b,theta=theta,T=T,D=D,R=R,nu=nu)
    ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    
    chi_range_c = np.arange(ulti_chi-np.radians(0.0001),ulti_chi+np.radians(0.0001),0.0000001)
    scores_range_c=magic_condition_angles(chi_range_c,theta=theta,T=T,D=D,R=R,nu=nu)
    ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    
    return ulti_chi

def mc_focus_solver(theta=None, chi = None, R=None, nu=None,D=None,T=None):
    
    chi_range_a = np.radians(np.arange(-80.0,80.0,0.01))
    scores_range_a = magic_condition_foci(chi_range_a,theta=theta,D=D,R=R,nu=nu)
    ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    plt.figure()
    plt.plot(chi_range_a,scores_range_a)
    plt.ylim(-10000000,10000000)
    plt.grid()
    
    chi_range_b = np.arange(ulti_chi-np.radians(1),ulti_chi+np.radians(1),0.0001)
    scores_range_b=magic_condition_foci(chi_range_b,theta=theta,D=D,R=R,nu=nu)
    ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    
    chi_range_c = np.arange(ulti_chi-np.radians(0.0001),ulti_chi+np.radians(0.0001),0.0000001)
    scores_range_c=magic_condition_foci(chi_range_c,theta=theta,D=D,R=R,nu=nu)
    ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    return ulti_chi

def radius_from_geo_focus(f_g, chi,theta, D):
    chi = np.radians(chi)
    theta=np.radians(theta)
    R = 2/(np.cos(chi-theta)/f_g - np.cos(chi+theta)/D)
    return R

def mc_solver_chi_focus(chi = None, f_g=None, theta=None, nu=None,D=None,T=None):
    ''' Given f_g (instead of R), solve the chi angle for magic condition
    '''
    chi_range_a = np.arange(-90.0,90.0,0.01)
    R = radius_from_geo_focus(f_g, chi_range_a, theta=theta, D=D)
    scores_range_a =  magic_condition_angles(chi_range_a,theta=theta,T=T,D=D,R=R,nu=nu)
    magic_chi = chi_range_a[abs(scores_range_a).argmin()]
#     plt.plot(chi_range_a,scores_range_a)
#     plt.ylim(-0.0005,0.0005)
#     plt.grid()
    
    chi_range_b = np.arange(magic_chi-1,magic_chi+1,0.0001)
    R = radius_from_geo_focus(f_g, chi_range_b, theta=theta, D=D)
    scores_range_b=magic_condition_angles(chi_range_b,theta=theta,T=T,D=D,R=R,nu=nu)
    magic_chi = chi_range_b[abs(scores_range_b).argmin()]
    
    chi_range_c = np.arange(magic_chi-0.0001,magic_chi+0.0001,0.0000001)
    R = radius_from_geo_focus(f_g, chi_range_c, theta=theta, D=D)
    scores_range_c=magic_condition_angles(chi_range_c,theta=theta,T=T,D=D,R=R,nu=nu)
    magic_chi = chi_range_c[abs(scores_range_c).argmin()]
    
    return magic_chi


def magic_condition_bend_radius(chi,theta,D,nu):
    R = ((1+nu)*D*np.sin(2*chi)*np.cos(chi-theta)-2*D*np.sin(theta-chi))/np.sin(2*theta)
    return R

def magic_condition_source_distance(chi,theta,R,nu):
    D = R*np.sin(2*theta)/((1+nu)*np.sin(2*chi)*np.cos(chi-theta)+2*np.sin(chi-theta))
    return D


def energy_resolution(chi,theta,R,D,T,nu):
    
    theta = np.radians(theta)
    chi = np.radians(chi)    
    term1 = T*(np.cos(chi)**2-nu*np.sin(chi)**2)/R
    term2 = T*(np.cos(theta)**2)/(2*D*np.cos(chi-theta))
    resolution = term1+term2
    return resolution




def ulti_magi_condi(chi,theta=8.99,nu=0.2,T=0.3,D=2200,return_all=False):
    theta = np.radians(theta)
    chi = np.radians(chi)
    R = ((1+nu)*D*np.sin(2*chi)*np.cos(chi-theta)-2*D*np.sin(theta-chi))/np.sin(2*theta)
    try: len(R)
    except: print('R:',R)
    # energy resolution equation
    term1 = T*(np.cos(chi)**2-nu*np.sin(chi)**2)/R
    term2 = T*(np.cos(theta)**2)/(2*D*np.cos(chi-theta))
    resolution = term1+term2
    if return_all:
        return resolution,term1,term2
    else: # so that we can use this function with `fsolve`
        return resolution

def ulti_magi_condi_with_R(chi,theta=8.99,nu=0.2,T=0.3,R=2000,return_all=False):
    theta = np.radians(theta)
    chi = np.radians(chi)
    D = R*np.sin(2*theta)/((1+nu)*np.sin(2*chi)*np.cos(chi-theta)+2*np.sin(chi-theta))
    try: len(D)
    except: print('D:',D)
    # energy resolution equation
    term1 = T*(np.cos(chi)**2-nu*np.sin(chi)**2)/R
    term2 = T*(np.cos(theta)**2)/(2*D*np.cos(chi-theta))
#     term3 = term1*term2
    resolution = term1+term2
    if return_all:
        return resolution,term1,term2
    else: # so that we can use this function with `fsolve`
        return resolution


def ulti_magi_condi_chi_solver(theta=8.98,T=0.3,D=22000,nu=0.2):
    chi_range_a = np.arange(-45.0,45.0,0.01)
    scores_range_a = ulti_magi_condi(chi_range_a,theta=theta,T=T,D=D,nu=nu)
    ulti_chi = chi_range_a[abs(scores_range_a).argmin()]
    plt.plot(chi_range_a,scores_range_a)
    plt.ylim(-0.0001,0.0001)
    plt.ylabel('$\Delta E/E$')
    plt.grid()
    chi_range_b = np.arange(ulti_chi-1,ulti_chi+1,0.0001)
    scores_range_b=ulti_magi_condi(chi_range_b,theta=theta,T=T,D=D,nu=nu)
    ulti_chi = chi_range_b[abs(scores_range_b).argmin()]
    
#     chi_range_c = np.arange(ulti_chi-0.0001,ulti_chi+0.0001,0.0000001)
#     scores_range_c=ulti_magi_condi(chi_range_c,theta=theta,T=T,D=D,nu=nu)
#     ulti_chi = chi_range_c[abs(scores_range_c).argmin()]
    
    return ulti_chi


def energy_divergence(theta1,theta2,hkl=[1,1,1]):
    import math_physics as mp
    energy1 = mp.bragg(hkl=hkl,theta=theta1)[0]
    energy2 = mp.bragg(hkl=hkl,theta=theta2)[0]
    energy_divergence = (energy2-energy1)*1000 # return in eV
    return energy_divergence


def dy_to_theta(dy,theta,chi,R,D=np.Infinity):
    theta=np.radians(theta)
    chi= np.radians(chi)
    y0 = -np.sin(theta+chi)*R
    theta_dy = np.arcsin((-y0-dy)/R)-chi-np.arcsin(dy/D)
    theta_dy = np.degrees(theta_dy)
    return theta_dy


def mc_y_evaluation(chi= 4.4671,theta=8.99,nu=0.2,T=0.3,R=2000,D=22000,d_y = np.arange(-2.5,2.5,0.001),verbose=False):
    theta_0=np.radians(theta)
    chi= np.radians(chi)
    y_0 = -np.sin(theta_0+chi)*R
    print(y_0)
    
    d_y = d_y
    theta = np.arcsin((-y_0-d_y)/R)-chi-np.arcsin(d_y/D) # radians
    mc_angle_distance = magic_condition_angle(np.degrees(chi),np.degrees(theta),nu,T,R,D) # in radians
    mbw = mono_beam_width(np.degrees(chi),np.degrees(theta),nu,T,R,D)
    if verbose:
        # Plot
        plt.figure(figsize=(13,5))
        plt.subplot(1,2,1)
        plt.plot(d_y,np.degrees(mc_angle_distance),label='Opening angle at Exit')
        plt.legend()
        plt.ylabel('deg')
        plt.xlabel('Relative vertical position (mm)')
        plt.title('Opening angle vs Vertical x-ray position')
        plt.subplot(1,2,2)
        plt.plot(d_y,mbw,color='r',label='Monochromatic beam width (mm)' )
        plt.xlabel('Relative vertical position (mm)')
        plt.legend()
        plt.title('Monochromatic beam width vs Vertical x-ray position')

    # Theta angle divergence over the whole beam vertically
    try:
        D_theta = np.degrees(theta[-1]-theta[0]) # in degree

        if verbose:
            print('Incident theta angle divergence over the whole beam vertically:\n',abs(D_theta))
            print('The biggest exiting point angle divergence is (deg):\n', np.degrees(max(abs(mc_angle_distance[-1]),abs(mc_angle_distance[0]))))
            D_theta_2 = 5/(R*np.cos(chi+theta_0))+5/D
            print(np.degrees(D_theta_2))
#             print(np.radians(D_theta)/np.tan(theta_0)*12.658)
        D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0]-mp.bragg(theta=np.degrees(theta[0]))[0]
        if verbose:
            print('Overall energy range:\n',D_energy)
    except:
        pass


    
    return mc_angle_distance, mbw


def mc_theta_evaluation(chi= 4.4671,theta=8.99,nu=0.2,T=0.3,R=2000,D=22000,verbose=False):
    theta_0=np.radians(theta)
    chi= np.radians(chi)

    d_theta = np.arange(-np.pi/4,np.pi/4,0.0001)
    theta = theta_0+d_theta # in radians
    mc_angle_distance = magic_condition_angle(np.degrees(chi),np.degrees(theta),nu,T,R,D) # in radians
    mbw = mono_beam_width(np.degrees(chi),np.degrees(theta),nu,T,R,D)
    if verbose:
        # Plot
#         plt.figure(figsize=(13,5))
#         plt.subplot(1,2,1)
#         plt.plot(np.degrees(d_theta),np.degrees(mc_angle_distance),label='Opening angle at Exit')
#         plt.legend()
#         plt.ylabel('deg')
#         plt.xlabel('Angle distance from Magic Condition Bragg angle (deg)')
#         plt.title('Opening angle vs Bragg angle')
#         plt.subplot(1,2,2)
        plt.plot(np.degrees(d_theta),mbw*1000,color='r',label='Monochromatic beam width ($\mu m$)' )
        plt.xlabel(r'Relative distance from magic condition Bragg angle $\theta_B$ (deg)')
        plt.ylabel('Monochromatic beam width ($\mu m$)')
#         plt.legend()
        plt.title('Monochromatic beam width vs Bragg angle')
        plt.grid()

#     # Theta angle divergence over the whole beam vertically
#     D_theta = np.degrees(d_theta[-1]-d_theta[0]) # in degree
#     if verbose:
#         print('Incident theta angle divergence over the whole beam vertically:\n',abs(D_theta))
#         print('The biggest exiting point angle divergence is (deg):\n', np.degrees(max(abs(mc_angle_distance[-1]),abs(mc_angle_distance[0]))))


#     D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0]-mp.bragg(theta=np.degrees(theta[0]))[0]
# #     print(theta[-1],theta[0])
#     if verbose:
#         print('Overall energy range: (keV)\n',D_energy)
    
    return mc_angle_distance, None, None


def mc_chi_evaluation(chi= 4.4671,theta=8.99,nu=0.2,T=0.3,R=2000,D=22000,verbose=False):
    theta=np.radians(theta)
    chi_0= np.radians(chi)

    d_chi = np.arange(-np.pi/4,np.pi/4,0.0001)
    chi = chi_0+d_chi # in radians
    mc_angle_distance = magic_condition_angle(np.degrees(chi),np.degrees(theta),nu,T,R,D) # in radians
    mbw = mono_beam_width(np.degrees(chi),np.degrees(theta),nu,T,R,D)
    if verbose:
#         # Plot
#         plt.figure(figsize=(13,5))
#         plt.subplot(1,2,1)
#         plt.plot(np.degrees(d_chi),np.degrees(mc_angle_distance),label='Opening angle at Exit')
#         plt.legend()
#         plt.ylabel('deg')
#         plt.xlabel('Angle distance from Magic Condition Chi angle (deg)')
#         plt.title('Opening angle VS asymmetry angle')
#         plt.subplot(1,2,2)
        plt.plot(np.degrees(d_chi),mbw*1000,color='r' )
        plt.xlabel('Relative distance from magic condition asymmetry angle $\chi$ (deg)')
        plt.ylabel('Monochromatic beam width ($\mu m$)')
#         plt.legend()
        plt.title('Monochromatic beam width VS asymmetry angle')
        plt.grid()

#     # Theta angle divergence over the whole beam vertically
#     D_theta = np.degrees(d_theta[-1]-d_theta[0]) # in degree
#     if verbose:
#         print('Incident theta angle divergence over the whole beam vertically:\n',abs(D_theta))
#         print('The biggest exiting point angle divergence is (deg):\n', np.degrees(max(abs(mc_angle_distance[-1]),abs(mc_angle_distance[0]))))


#     D_energy = mp.bragg(theta=np.degrees(theta[-1]))[0]-mp.bragg(theta=np.degrees(theta[0]))[0]
#     print(theta[-1],theta[0])
#     if verbose:
#         print('Overall energy range: (keV)\n',D_energy)
    
    return mc_angle_distance, None, None

def mc_theta_wanderoff(chi= 4.4671,theta=8.99,nu=0.2,T=0.3,R=2000,D=22000,verbose=False):
    theta_0=np.radians(theta)
    chi = np.radians(chi)
    d_theta = np.arange(-np.pi/4,np.pi/4,0.0001)
    theta = theta_0+d_theta # in radians
    mc_angle_distance = magic_condition_angles(np.degrees(chi),np.degrees(theta),nu,T,R,D) # in radians
    mbw = mono_beam_width(np.degrees(chi),np.degrees(theta),nu,T,R,D)
    return np.degrees(d_theta), mbw*1000
def mc_chi_wanderoff(chi= 4.4671,theta=8.99,nu=0.2,T=0.3,R=2000,D=22000,verbose=False):
    theta=np.radians(theta)
    chi_0= np.radians(chi)
    d_chi = np.arange(-np.pi/4,np.pi/4,0.0001)
    chi = chi_0+d_chi # in radians
    mc_angle_distance = magic_condition_angles(np.degrees(chi),np.degrees(theta),nu,T,R,D) # in radians
    mbw = mono_beam_width(np.degrees(chi),np.degrees(theta),nu,T,R,D)
    return np.degrees(d_chi), mbw*1000
