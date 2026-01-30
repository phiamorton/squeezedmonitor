import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


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
# Parameters
# ----------------------------
d = 1.0
y = np.linspace(-d, d, 1000)

# Non-focused beam: weak modulation
Ny_nonfocused = 1.0 + 0.08 * np.cos(2 * np.pi * y / d)

# Focused beam: strong modulation
Ny_focused = 1.0 + 0.9 * np.cos(2 * np.pi * y / d)

# ----------------------------
# Figure layout
# ----------------------------
fig = plt.figure(figsize=(10, 6))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 3, 0.2], height_ratios=[1, 1])

# ----------------------------
# Beam schematic (top)
# ----------------------------
ax_beam1 = fig.add_subplot(gs[0, 0])
for yy in np.linspace(-0.4, 0.4, 6):
    ax_beam1.plot([-1, 1], [yy, yy], color="gray", lw=5)
#ellipse1 = plt.Circle((0, 0), 0.25, color="black")
#ax_beam1.add_patch(ellipse1)
ax_beam1.add_patch(Ellipse((0, 0), width=1.5, height=.5,
                     edgecolor='black',
                     facecolor='black',
                     linewidth=1,zorder=5))
ax_beam1.set_xlim(-1, 1)
ax_beam1.set_ylim(-0.6, 0.6)
ax_beam1.axis("off")

# ----------------------------
# Beam schematic (bottom)
# ----------------------------
ax_beam2 = fig.add_subplot(gs[1, 0])
for yy in np.linspace(-0.4, 0.4, 6):
    ax_beam2.plot([-1, 1], [yy, yy], color="gray", lw=5)
#ellipse2 = plt.Circle((0, 0), 0.08, color="black")
#ax_beam2.add_patch(ellipse2)
ax_beam2.add_patch(Ellipse((0, 0), width=.32, height=.12,
                     edgecolor='black',
                     facecolor='black',
                     linewidth=1,zorder=5))
ax_beam2.annotate("", xy=(0.5, -0.1), xytext=(0.5, 0.1),
                  arrowprops=dict(arrowstyle="<->"))
ax_beam2.text(0.55, 0, "d", va="center")
ax_beam2.set_xlim(-1, 1)
ax_beam2.set_ylim(-0.6, 0.6)
ax_beam2.axis("off")

# ----------------------------
# Non-focused beam plot
# ----------------------------
ax1 = fig.add_subplot(gs[0, 1])
ax1.plot(y, Ny_nonfocused, color="black")
ax1.set_ylim(0, 2)
ax1.set_xlim(-d, d)
ax1.set_ylabel(r"$N_\gamma$ normalized")
ax1.set_xticks([-d/2, 0, d/2])
ax1.set_xticklabels([r"$-d/2$", "0", r"$d/2$"])
ax1.set_title("(a) Non focused beam", y=-0.4)
ax1.grid(True, alpha=0.3)

# ----------------------------
# Focused beam plot
# ----------------------------
ax2 = fig.add_subplot(gs[1, 1])
ax2.plot(y, Ny_focused, color="black")
ax2.set_ylim(0, 2)
ax2.set_xlim(-d, d)
ax2.set_ylabel(r"$N_\gamma$ normalized")
#ax2.set_xlabel("y")
ax2.set_xticks([-d/2, 0, d/2])
ax2.set_xticklabels([r"$-d/2$", "0", r"$d/2$"])
ax2.set_title("(b) Focused beam", y=-0.4) 
ax2.text(-0.05, 1.65, r"$N_{\max}$")
ax2.text(0.45*d, 0.25, r"$N_{\min}$")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
