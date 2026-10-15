"""Graphical abstract: two stacked reservoirs across the freezing front. Greyscale."""
import warnings
warnings.filterwarnings("ignore")
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle, FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
import json
_B = json.load(open(os.path.join(HERE, "production", "best_fit.json")))
_L, _S = _B["liquid_wt_pct"], _B["solid_wt_pct"]
LIQ_LABEL = f"O {_L['O']:.1f}, Si {_L['Si']:.1f}, S {_L['S']:.1f}, C {_L['C']:.1f} wt%"
SOL_LABEL = f"Fe Ni, C {_S['C']:.1f}, Si {_S['Si']:.1f}, S {_S['S']:.1f}, O {_S['O']:.2f} wt%"
os.makedirs(FIGDIR, exist_ok=True)

SOURCE = ("Source: ICB partition model with nickel (icb_partition_model); "
          "PREM Dziewonski & Anderson 1981; partition D_i Alfe 2002, Li 2019, "
          "Badro 2014; shear Murphy 2013, Ono 2013, Prescher 2015. "
          "Greyscale per house figure standard.")

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
})

# Elsevier graphical abstract ~ 13 x 7 cm, ratio ~1.857
fig = plt.figure(figsize=(13 / 2.54, 7 / 2.54))
ax = fig.add_axes([0.0, 0.0, 1.0, 1.0])
ax.set_xlim(0, 1.857)
ax.set_ylim(0, 1.0)
ax.axis("off")

W = 1.857

# Outer box (four spine thin black frame around canvas)
ax.add_patch(Rectangle((0.012, 0.012), W - 0.024, 1.0 - 0.024, fill=False,
                        edgecolor="black", linewidth=0.9))

# Reservoir geometry
left, right = 0.18, 1.857 - 0.18
top_band = (0.55, 0.86)     # liquid outer core (top)
bot_band = (0.20, 0.50)     # solid inner core (bottom)

# Labels sit on white boxes so the hatching cannot run through the text.
# Top reservoir: liquid outer core, mid grey, diagonal hatch
ax.add_patch(Rectangle((left, top_band[0]), right - left,
                       top_band[1] - top_band[0], facecolor="0.6",
                       edgecolor="black", linewidth=0.9, hatch="////"))
ax.text((left + right) / 2, (top_band[0] + top_band[1]) / 2 + 0.04,
        "OUTER CORE LIQUID", ha="center", va="center", fontweight="bold",
        fontsize=8.5, color="black", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
ax.text((left + right) / 2, (top_band[0] + top_band[1]) / 2 - 0.05,
        LIQ_LABEL, ha="center", va="center",
        fontsize=7.5, color="black", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

# Bottom reservoir: Fe Ni solid with carbon, silicon and sulfur, light grey, dotted hatch
ax.add_patch(Rectangle((left, bot_band[0]), right - left,
                       bot_band[1] - bot_band[0], facecolor="0.88",
                       edgecolor="black", linewidth=0.9, hatch="...."))
ax.text((left + right) / 2, (bot_band[0] + bot_band[1]) / 2 + 0.04,
        "INNER CORE SOLID", ha="center", va="center", fontweight="bold",
        fontsize=8.5, color="black", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))
ax.text((left + right) / 2, (bot_band[0] + bot_band[1]) / 2 - 0.05,
        SOL_LABEL, ha="center", va="center",
        fontsize=7.5, color="black", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

# Freezing front line between the two reservoirs
front_y = 0.525
ax.plot([left, right], [front_y, front_y], color="black", lw=1.3, ls="-")
ax.text(right + 0.01, front_y, "freezing\nfront", ha="left", va="center",
        fontsize=6.8, color="black")

# Partition arrows downward through the front, with D_i labels
import icb_partition_model as M
_D = {el: M.PP["partition_coefficients"][el]["D"] for el in ("O", "Si", "S", "C")}
for xc, lab in [(0.30, f"O  D {_D['O']:g}"), (0.62, f"Si  D {_D['Si']:g}"), (0.94, f"S  D {_D['S']:g}"),
                (1.26, f"C  D {_D['C']:g}")]:
    ar = FancyArrowPatch((xc, top_band[0] - 0.01), (xc, bot_band[1] + 0.01),
                         arrowstyle="-|>", mutation_scale=12, color="black",
                         lw=1.2, shrinkA=0, shrinkB=0)
    ax.add_patch(ar)
    ax.text(xc + 0.015, front_y, lab, ha="left", va="center", fontsize=6.6,
            color="black", bbox=dict(facecolor="white", edgecolor="none", pad=1.5))

out = os.path.join(FIGDIR, "graphical_abstract.png")
fig.savefig(out, dpi=600, metadata={"Software": None})
plt.close(fig)
print("wrote", out, "ratio", round(W, 3))
