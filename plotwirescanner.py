import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl

mpl.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 16,
    "axes.labelsize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
})
# ----------------------------
# Parameters (nm)
# ----------------------------
sigma_e = 100.0
sigma_L = 50.0
lambda0 = 800.0
wiggle_amp = sigma_L
y0 = 150.0

# ----------------------------
# Zoomed coordinate grids
# ----------------------------
x = np.linspace(-800, 800, 3000)
y = np.linspace(-800, 800, 2000)
X, Y = np.meshgrid(x, y, indexing="ij")

# ----------------------------
# Electron beam (2D Gaussian)
# ----------------------------
E_e = np.exp(-(X**2 + Y**2) / (2 * sigma_e**2))

# ----------------------------
# Laser envelope (opaque bar)
# ----------------------------
laser_env = np.exp(-(Y - y0)**2 / (2 * sigma_L**2))

# Single laser field line
k = 2 * np.pi / lambda0
laser_y = y0 + wiggle_amp * np.cos(k * x)

# ----------------------------
# Reconstructed signal (tilted)
# ----------------------------
signal_y = np.exp(-y**2 / (2 * (sigma_e**2 + sigma_L**2)))

# ----------------------------
# Plot
# ----------------------------
fig, axes = plt.subplots(
    1, 2, figsize=(9.5, 5),
    sharey=True,
    gridspec_kw={"width_ratios": [1.3, 0.7]}
)

# ==========================================================
# LEFT: Spatial schematic
# ==========================================================
ax = axes[0]

# Electron beam
ax.imshow(
    E_e.T,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin="lower",
    cmap="Blues",
    alpha=0.9,
    aspect="auto"
)

# Electron beam 1σ outline
ax.contour(
    X, Y, E_e,
    levels=[np.exp(-0.5)],
    colors="blue",
    linewidths=2
)

# Laser envelope as opaque bar
ax.imshow(
    laser_env.T,
    extent=[x.min(), x.max(), y.min(), y.max()],
    origin="lower",
    cmap="Reds",
    alpha=0.35,
    aspect="auto"
)

# Single laser field line
#ax.plot(x, laser_y, color="darkred", linewidth=2,label=f"{lambda0} nm laser focused to {sigma_L} nm")

ax.set_xlim(-800, 800)
ax.set_ylim(-800, 800)

ax.set_xlabel("x (nm)")
ax.set_ylabel("y (nm)")
ax.set_title("Laser–electron overlap (single position)")
ax.text(-150, -200, "Electron beam", color="blue", fontsize=14)
ax.text(-400, 300, "Laser beam", color="red", fontsize=14)
# ==========================================================
# RIGHT: Tilted reconstructed profile
# ==========================================================
ax2 = axes[1]

ax2.plot(signal_y, y, color="black", linewidth=2)
ax2.fill_betweenx(y, 0, signal_y, color="gray", alpha=0.3)
ax2.axhline(y0, color="red", linestyle="--", alpha=0.7)
ax2.set_xlabel("Compton signal (a.u.)")
ax2.set_title("Reconstructed profile")
ax2.grid(True, alpha=0.3)

# ----------------------------
# Force identical y-scale
# ----------------------------
ax.set_ylim(y.min(), y.max())
ax2.set_ylim(y.min(), y.max())

plt.tight_layout()
plt.show()

