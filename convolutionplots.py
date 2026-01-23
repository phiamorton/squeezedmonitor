import numpy as np
import matplotlib.pyplot as plt
import cmath

# ----------------------------
# Parameters
# ----------------------------
phase_dep=True
dB_squeezing=10
A = 1.0
lambda_nm = 800.0      # nm

num_to_avg = 10
dA_noise = 0.05*A
dphi_noise = dA_noise

dx_beam = 100.0  # nm
N_beam = 100.0
     
# coordinate
x = np.linspace(-lambda_nm, lambda_nm, 2000)
dx = x[1] - x[0]

# ----------------------------
# Laser fringe
# ----------------------------
def laser_fringe(x, A, lambd, phi):
    return A**2 * np.sin(2*np.pi*x /lambd + phi)**2

# ----------------------------
# Gaussian beam
# ----------------------------
def gaussian_beam(x, N, dx):
    return N  * np.exp(-0.5 * (x/dx)**2) *1/ (np.sqrt(2*np.pi)*dx)

#beam = gaussian_beam(x, N_beam, dx_beam)
#plt.plot(x,laser_fringe(x, A, lambda_nm, 0))
#plt.show()

# ============================================================
# Phase scan with averaging
# ============================================================

num_phase = 1000
num_avg   = 50
phi_scan = np.linspace(0, 6.28, num_phase)

mean_both = np.zeros_like(phi_scan)
std_both  = np.zeros_like(phi_scan)

mean_A = np.zeros_like(phi_scan)
std_A  = np.zeros_like(phi_scan)

def phiscan(beamsize):

    beam = gaussian_beam(x, N_beam, beamsize)
    mean_both = np.zeros_like(phi_scan)
    std_both  = np.zeros_like(phi_scan)
    mean_A    = np.zeros_like(phi_scan)
    std_A     = np.zeros_like(phi_scan)
    mean_p    = np.zeros_like(phi_scan)
    std_p     = np.zeros_like(phi_scan)

    for i, phi in enumerate(phi_scan):
        signals_both = []
        signals_A    = []
        signals_p    = []
        ###HERE NEED TO FIGURE OUT SQUEEZING ANGLE
        if phase_dep==True:
            squeezing_parameter= 10**(-dB_squeezing*(2*np.sin(2*np.pi*(phi)/lambda_nm+phi)**2-1)/20)+0.00001
            # here is the amount of squeezing and phase for squeezing quadrature 
        else:
            squeezing_parameter= 1/2*np.log(10**(-dB_squeezing/10))
            r=np.log(dB_squeezing)/2
            #print(squeezing_parameter)
            squeezing_parameter=10**(-dB_squeezing/20)
        for _ in range(num_avg):
            laser_fringe_weight=laser_fringe(x,A,lambda_nm, phi)
            dA   = np.random.normal(0, dA_noise)
            dphi = np.random.normal(0, dphi_noise)
            #print(dA,dphi)
            #print(dA,dphi)
            
            fringe = laser_fringe(x, A + dA, lambda_nm, phi+dphi)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_both.append(np.trapezoid(signal_x, x))
            #signals_both.append(fringe)
            #print((A+dA)*(phi+dphi))

            fringe = laser_fringe(x, A + dA/squeezing_parameter, lambda_nm, phi+dphi*squeezing_parameter)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_A.append(np.trapezoid(signal_x, x))
            #signals_A.append(fringe)
            #print((A+dA/squeezing_parameter)*(phi+dphi*squeezing_parameter))

            fringe = laser_fringe(x, A + dA*squeezing_parameter, lambda_nm, phi+ dphi/squeezing_parameter)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_p.append(np.trapezoid(signal_x, x))
            #print((A+dA*squeezing_parameter)*(phi+dphi/squeezing_parameter))
            #signals_p.append(fringe)

        mean_both[i] = np.mean(signals_both)
        std_both[i]  = np.std(signals_both)
        mean_A[i]    = np.mean(signals_A)
        std_A[i]     = np.std(signals_A)
        mean_p[i]    = np.mean(signals_p)
        std_p[i]     = np.std(signals_p)

    total_variance = np.trapezoid(std_both, phi_scan)
    total_variance_A = np.trapezoid(std_A, phi_scan)
    total_variance_p = np.trapezoid(std_p, phi_scan)
    print(total_variance, total_variance_A, total_variance_p)

    return mean_both, std_both, mean_A, std_A, mean_p, std_p

def plotphiscan(beamsize): 
    mean_both1, std_both1, mean_A1, std_A1, mean_p, std_p=phiscan(beamsize)

    if phase_dep==True:
        plt.plot(phi_scan, mean_p, lw=2, label=f"with {dB_squeezing} dB relative phase dependent squeezing, {beamsize} nm e beam")
        plt.fill_between(
            phi_scan, mean_p - std_p, mean_p + std_p, alpha=0.3
        )
    else:
        """
        plt.plot(phi_scan, mean_A1, lw=2, label=f"with {dB_squeezing} dB phase squeezing, {beamsize} nm e beam")
        plt.fill_between(
            phi_scan, mean_A1 - std_A1, mean_A1 + std_A1, alpha=0.3
        )
        """
        plt.plot(phi_scan, mean_p, lw=2, label=f"with {dB_squeezing} dB amplitude squeezing, {beamsize} nm e beam")
        plt.fill_between(
            phi_scan, mean_p - std_p, mean_p + std_p, alpha=0.3
        )
        #print(std_both1[500], std_A1[500])
    
    plt.plot(phi_scan, mean_both1, lw=2, label=f"without squeezing, {beamsize} nm e beam")
    plt.fill_between(
        phi_scan, mean_both1 - std_both1, mean_both1 + std_both1, alpha=0.3
    )
    
    
    

def ideal_fringe(beamsize):
    beam = gaussian_beam(x, N_beam, beamsize)
    signal = []
    for phi in enumerate(phi_scan):
        ideal_fringe=laser_fringe(x,A,lambda_nm, phi[1])
        signal_x = ideal_fringe*beam 
        signal.append(np.trapezoid(signal_x, x))
    return signal 
    

plt.figure(figsize=(7,4))
beamsize1=10
plotphiscan(beamsize1)
signal1=ideal_fringe(beamsize1)
beamsize2=50
plotphiscan(beamsize2)
signal2=ideal_fringe(beamsize2)
#plt.plot(phi_scan,signal1, label=f'fringe pattern with 0 error, {beamsize1} nm beam', color='brown')
#plt.plot(phi_scan,signal2, label=f'fringe pattern with 0 error, {beamsize2} nm beam', color='black')
plt.xlabel("Relative phase between beam and fringe $\\phi$")
plt.ylabel("Detected signal (a.u.)")
plt.title("Phase scan: amplitude noise vs phase + amplitude noise")
plt.legend()
plt.tight_layout()
#plt.xlim([np.pi-np.pi/8,np.pi+np.pi/8])
#plt.ylim([0,4000])
plt.title(f'Detected photons averaged over {num_avg} shots with {lambda_nm} nm fringe spacing and {dA_noise*100}% laser noise')
plt.show()

#plot difference
"""
mean_both1, std_both1, mean_A1, std_A1, mean_p1, std_p1=phiscan(beamsize1)
mean_both2, std_both2, mean_A2, std_A2,mean_p2, std_p2=phiscan(beamsize2)
plt.figure(figsize=(7,4))

plt.plot(phi_scan, np.abs(mean_p1-mean_p2), label='amplitude squeezed')
plt.fill_between(
    phi_scan, np.abs(mean_p1-mean_p2) - np.sqrt(std_p1**2+std_p2**2),np.abs(mean_p1-mean_p2) + np.sqrt(std_p1**2+std_p2**2), alpha=0.3
)
plt.plot(phi_scan, np.abs(mean_A1-mean_A2), label='phase squeezed')
plt.fill_between(
    phi_scan, np.abs(mean_A1-mean_A2) - np.sqrt(std_A1**2+std_A2**2),np.abs(mean_A1-mean_A2) + np.sqrt(std_A1**2+std_A2**2), alpha=0.3
)
plt.plot(phi_scan, np.abs(mean_both1-mean_both2), label='not squeezed')
plt.fill_between(
    phi_scan, np.abs(mean_both1-mean_both2) - np.sqrt(std_both1**2+std_both2**2),np.abs(mean_both1-mean_both2) + np.sqrt(std_both1**2+std_both2**2), alpha=0.3
)
plt.xlabel("Relative phase $\\phi$")
plt.ylabel(f"Detected signal difference between {beamsize1} nm and {beamsize2} nm beam (arb.)")
plt.xlim([np.pi-np.pi/8,np.pi+np.pi/8])
plt.title(f'Detected photons difference averaged over {num_avg} shots')
plt.legend()
plt.show()
"""



