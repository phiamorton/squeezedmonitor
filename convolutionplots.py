import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Parameters
# ----------------------------
A = 1.0
lambda_nm = 800.0      # nm

num_to_avg = 10
dA_noise = 0.01
dphi_noise = 0.01

dx_beam = 100.0  # nm
N_beam = 1000.0

# coordinate
x = np.linspace(-5000, 5000, 20001)
dx = x[1] - x[0]

# ----------------------------
# Laser fringe
# ----------------------------
def laser_fringe(x, A, lambd, phi):
    return A**2 * np.sin(x / lambd / 2 + phi)**2

# ----------------------------
# Gaussian beam
# ----------------------------
def gaussian_beam(x, N, dx):
    return N  * np.exp(-0.5 * (x/dx)**2) *1/ (np.sqrt(2*np.pi)*dx)

beam = gaussian_beam(x, N_beam, dx_beam)

# ============================================================
# Spatial averages: (dA + dphi) vs (dA only)
# ============================================================
"""
phi0 = 0.0

fringe_avg_both = np.zeros_like(x)
fringe_avg_A    = np.zeros_like(x)

for _ in range(num_to_avg):
    # both noises
    dA   = np.random.uniform(-dA_noise, dA_noise)
    dphi = np.random.uniform(-dphi_noise, dphi_noise)
    fringe_avg_both += laser_fringe(x, A + dA, lambda_nm, phi0 + dphi)

    # amplitude only
    dA = np.random.uniform(-dA_noise, dA_noise)
    fringe_avg_A += laser_fringe(x, A + np.sqrt(dA**2+dphi**2), lambda_nm, phi0)

fringe_avg_both /= num_to_avg
fringe_avg_A    /= num_to_avg

signal_both = np.convolve(fringe_avg_both, beam,mode="same")
signal_A    = np.convolve(fringe_avg_A, beam,mode="same")

# ----------------------------
# Plot spatial comparison
# ----------------------------
plt.figure(figsize=(8,5))
#plt.plot(x, fringe_avg_both, label="⟨Fringe⟩ (dA + dφ)")
#plt.plot(x, fringe_avg_A, '--', label="⟨Fringe⟩")
#plt.plot(x, signal_both, label="Not Squeezed signal", lw=2)
#plt.plot(x, signal_A, '--', label="Squeezed Signal", lw=2)
#plt.plot(x, gaussian_beam(x, N_beam, 100), label="100 nm beam")
#plt.plot(x, gaussian_beam(x, N_beam, 1), label="1 nm beam")
size1=1
size2=100
size3=10
plt.plot(x, np.convolve(fringe_avg_A, gaussian_beam(x, N_beam, size1),mode="same"), label=f"{size1} nm beam, squeezed")
plt.plot(x, np.convolve(fringe_avg_both, gaussian_beam(x, N_beam, size1),mode="same"), label=f"{size1} nm beam, not squeezed")
plt.plot(x, np.convolve(fringe_avg_A, gaussian_beam(x, N_beam, size2),mode="same"), label=f"{size2} nm beam, squeezed")
plt.plot(x, np.convolve(fringe_avg_both, gaussian_beam(x, N_beam, size2),mode="same"), label=f"{size2} nm beam, not squeezed")
plt.plot(x, np.convolve(fringe_avg_A, gaussian_beam(x, N_beam, size3),mode="same"), label=f"{size3} nm beam, squeezed")
plt.plot(x, np.convolve(fringe_avg_both, gaussian_beam(x, N_beam, size3),mode="same"), label=f"{size3} nm beam, not squeezed")
plt.xlabel("x [nm]")
plt.ylabel("Intensity of Compton Scattered Photons (arb.)")
plt.legend()
plt.title(f'Detected photons at the fringe minima averaged over {num_to_avg} shots')
plt.xlim([-lambda_nm/8,lambda_nm/8])
plt.ylim([0,15])
plt.tight_layout()
#plt.show()
"""
# ============================================================
# Phase scan with averaging
# ============================================================

num_phase = 1000
num_avg   = 500
squeezing_parameter=0.3
#phi_scan = np.linspace(0, 2*np.pi, num_phase)
phi_scan = np.linspace(3, 3.3, num_phase)

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

        for _ in range(num_avg):
            laser_fringe_weight=laser_fringe(x,A,lambda_nm, phi)
            dA   = np.random.normal(0, dA_noise)
            dphi = np.random.normal(0, dphi_noise)
            #print(dA,dphi)
            
            fringe = laser_fringe(x, A + dA, lambda_nm, phi + dphi)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_both.append(np.trapezoid(signal_x, x))
            #signals_both.append(fringe)

            fringe = laser_fringe(x, A + dA/squeezing_parameter, lambda_nm, phi+dphi*squeezing_parameter)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_A.append(np.trapezoid(signal_x, x))
            #signals_A.append(fringe)

            fringe = laser_fringe(x, A + dA*squeezing_parameter, lambda_nm, phi+ dphi/squeezing_parameter)
            signal_x = fringe*beam #np.convolve(fringe, beam, mode="same") * dx
            signals_p.append(np.trapezoid(signal_x, x))
            #signals_p.append(fringe)
            

        mean_both[i] = np.mean(signals_both)
        std_both[i]  = np.std(signals_both)
        mean_A[i]    = np.mean(signals_A)
        std_A[i]     = np.std(signals_A)
        mean_p[i]    = np.mean(signals_p)
        std_p[i]     = np.std(signals_p)

    total_variance = np.trapezoid(std_both**2, phi_scan)
    total_variance_A = np.trapezoid(std_A**2, phi_scan)
    total_variance_p = np.trapezoid(std_p**2, phi_scan)
    print(total_variance, total_variance_A, total_variance_p)

    return mean_both, std_both, mean_A, std_A, mean_p, std_p

def plotphiscan(beamsize):
    mean_both1, std_both1, mean_A1, std_A1, mean_p, std_p=phiscan(beamsize)
    
    plt.plot(phi_scan, mean_A1, lw=2, label=f"with phase squeezing, {beamsize} nm beam")
    plt.fill_between(
        phi_scan, mean_A1 - std_A1, mean_A1 + std_A1, alpha=0.3
    )

    plt.plot(phi_scan, mean_p, lw=2, label=f"with amplitude squeezing, {beamsize} nm beam")
    plt.fill_between(
        phi_scan, mean_p - std_p, mean_p + std_p, alpha=0.3
    )
    #print(std_both1[500], std_A1[500])
    
    plt.plot(phi_scan, mean_both1, lw=2, label=f"without squeezing, {beamsize} nm beam")
    plt.fill_between(
        phi_scan, mean_both1 - std_both1, mean_both1 + std_both1, alpha=0.3
    )

plt.figure(figsize=(7,4))
beamsize1=100
plotphiscan(beamsize1)
#plotphiscan(50)
beamsize2=10
#plotphiscan(beamsize2)
plt.xlabel("Relative phase between beam and fringe $\\phi$")
plt.ylabel("Detected signal (a.u.)")
plt.title("Phase scan: amplitude noise vs phase + amplitude noise")
plt.legend()
plt.tight_layout()
#plt.xlim([np.pi-np.pi/8,np.pi+np.pi/8])
#plt.ylim([0,15])
plt.title(f'Detected photons averaged over {num_avg} shots')
plt.show()

#difference
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



