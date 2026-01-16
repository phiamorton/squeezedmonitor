import numpy as np
import matplotlib.pyplot as plt

# ----------------------------
# Parameters
# ----------------------------
A = 1.0
lambda_nm = 800.0      # nm
phi = 0.0

num_to_avg = 1000
dA_noise = 0.1
dphi_noise = 0

dx_beam = 10.0 # nm
N_beam = 1.0

# coordinate
x = np.linspace(-5000, 5000, 20001)  # nm
dx = x[1] - x[0]

# ----------------------------
# Laser fringe
# ----------------------------
def laser_fringe(x, A, lambd, phi):
    return A**2 * np.sin(x / lambd / 2 + phi)**2

# ----------------------------
# Average over noise
# ----------------------------
fringe_avg = np.zeros_like(x)

for _ in range(num_to_avg):
    dA = np.random.uniform(-dA_noise, dA_noise)
    dphi = np.random.uniform(-dphi_noise, dphi_noise)
    fringe_avg += laser_fringe(x, A + dA, lambda_nm, phi + dphi)

fringe_avg /= num_to_avg

# ----------------------------
# Gaussian beam
# ----------------------------
def gaussian_beam(x, N, dx):
    return N / (np.sqrt(2*np.pi)*dx) * np.exp(-0.5 * (x/dx)**2)

beam = gaussian_beam(x, N_beam, dx_beam)

# ----------------------------
# Convolution
# ----------------------------
signal = np.convolve(fringe_avg, beam, mode="same") * dx

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(8,5))
plt.plot(x, fringe_avg, label="Averaged fringe")
plt.plot(x, beam / beam.max() * fringe_avg.max(), '--', label="Beam (scaled)")
plt.plot(x, signal, label="Convolved signal", linewidth=2)
plt.xlabel("x [nm]")
plt.ylabel("Intensity (arb.)")
plt.legend()
plt.tight_layout()
plt.show()

# ----------------------------
# Phase scan with averaging
# ----------------------------

num_phase = 50       # phase scan points
num_avg = 50         # averages per phase
phi_scan = np.linspace(0, 2*np.pi, num_phase)

signal_mean = np.zeros_like(phi_scan)
signal_std  = np.zeros_like(phi_scan)

for i, phi in enumerate(phi_scan):

    signals = []

    for _ in range(num_avg):
        dA   = np.random.uniform(-dA_noise, dA_noise)
        dphi = np.random.uniform(-dphi_noise, dphi_noise)

        fringe = laser_fringe(x, A + dA, lambda_nm, phi + dphi)
        signal_x = np.convolve(fringe, beam, mode="same") * dx

        # detector integrates over x
        signals.append(np.trapz(signal_x, x))

    signals = np.array(signals)
    signal_mean[i] = signals.mean()
    signal_std[i]  = signals.std()

# ----------------------------
# Plot
# ----------------------------
plt.figure(figsize=(7,4))
plt.plot(phi_scan, signal_mean, lw=2, label="Mean signal")
plt.fill_between(
    phi_scan,
    signal_mean - signal_std,
    signal_mean + signal_std,
    alpha=0.3,
    label=r"$\pm1\sigma$"
)

plt.xlabel("Relative phase $\\phi$")
plt.ylabel("Detected signal (arb.)")
plt.title("Phase scan with averaging at each phase point")
plt.legend()
plt.tight_layout()
plt.show()
