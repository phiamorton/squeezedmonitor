import numpy as np
import matplotlib.pyplot as plt

plt.figure(figsize=(8,5))
beams=[10,20]
grid_width=np.max(beams)*3
for beamsize in beams:
    # ----------------------------
    # Beam parameters (nm)
    # ----------------------------
    sigma_e = beamsize     # electron beam rms
    sigma_L = 50.0      # laser rms
    lambda0=800 #nm laser wavelength
    # ----------------------------
    # Noise parameters
    # ----------------------------
    A0 = 1.0
    sigma_A0 = 0.05
    sigma_phi0 = 0.05
    r = 1.0              # squeezing factor

    # ----------------------------
    # Coordinate grid (nm)
    # ----------------------------
    x = np.linspace(-grid_width, grid_width, 2001)
    y = np.linspace(-grid_width, grid_width, 2001)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    X, Y = np.meshgrid(x, y, indexing="ij")

    # ----------------------------
    # Scan & sampling
    # ----------------------------
    y_scan = np.linspace(-grid_width, grid_width, 150)
    N_samples = 50

    # ----------------------------
    # Electron beam
    # ----------------------------
    E_e = 1/(np.sqrt(2*np.pi)*sigma_e**2)*np.exp(-(X**2 + Y**2) / (2 * sigma_e**2))

    # ----------------------------
    # Overlap
    # ----------------------------
    def overlap(Ee, EL):
        return np.sum(Ee * EL) * dx * dy

    # ----------------------------
    # Noise cases
    # ----------------------------
    cases = {
        "Unsqueezed": (sigma_A0, sigma_phi0),
        "Amplitude-squeezed": (sigma_A0*np.exp(-r), sigma_phi0*np.exp(+r)),
        "Simulated beam": (0,0)
    }

    results = {}

    # ----------------------------
    # Monte Carlo scan
    # ----------------------------
    for label, (sigma_A, sigma_phi) in cases.items():

        mean_vals = []
        std_vals = []

        for y0 in y_scan:

            vals = []

            for _ in range(N_samples):

                A = A0 * (1 + np.random.normal(0, sigma_A))
                phi = np.random.normal(0, sigma_phi)
                k = 2*np.pi/lambda0
                E_L = (
                    A * np.cos(k*X+phi)
                    * np.exp(-((Y - y0)**2) / (2 * sigma_L**2))
                )

                vals.append(overlap(E_e, E_L))

            mean_vals.append(np.mean(vals))
            std_vals.append(np.std(vals))

        mean_vals = np.array(mean_vals)
        std_vals = np.array(std_vals)

        # Normalize to compare profiles
        #mean_vals /= mean_vals.max()
        #std_vals /= mean_vals.max()

        results[label] = (mean_vals, std_vals)

    # ----------------------------
    # Plot
    # ----------------------------

    for label, (mean_vals, std_vals) in results.items():
        if label=="Simulated beam":
            label=f"simulated beam, rms size {sigma_e} nm"
        plt.plot(y_scan, mean_vals, label=label)
        plt.fill_between(
            y_scan,
            mean_vals - std_vals,
            mean_vals + std_vals,
            alpha=0.3
        )
#plt.plot(x,np.exp(-(x**2) / (2 * sigma_e**2)),label=f'simulated electron beam')#, $\sigma_x$={sigma_e} nm')
plt.xlabel("laser position $y_0$  (nm)")
plt.ylabel("compton signal [a.u.]")
plt.title(f"reconstruction of an electron beam transverse profile from {N_samples} samples\nwith a {lambda0} nm laser & {sigma_A0*100} % noise focused to {sigma_L} nm rms")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
