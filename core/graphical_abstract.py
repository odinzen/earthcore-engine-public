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

# Title
ax.text(W / 2, 0.93, "OXYGEN RICH OUTER CORE FEEDS A NEARLY PURE Fe Ni SOLID",
        ha="center", va="center", fontweight="bold", fontsize=9.5,
        color="black")

# Reservoir geometry
left, right = 0.18, 1.857 - 0.18
top_band = (0.55, 0.86)     # liquid outer core (top)
bot_band = (0.20, 0.50)     # solid inner core (bottom)

# Top reservoir: oxygen rich liquid outer core, mid grey, diagonal hatch
ax.add_patch(Rectangle((left, top_band[0]), right - left,
                       top_band[1] - top_band[0], facecolor="0.6",
                       edgecolor="black", linewidth=0.9, hatch="////"))
ax.text((left + right) / 2, (top_band[0] + top_band[1]) / 2 + 0.04,
        "OUTER CORE LIQUID", ha="center", va="center", fontweight="bold",
        fontsize=8.5, color="black")
ax.text((left + right) / 2, (top_band[0] + top_band[1]) / 2 - 0.05,
        "O 18, S 1.5 mol %, plus Ni 5", ha="center", va="center",
        fontsize=7.5, color="black")

# Bottom reservoir: nearly pure Fe Ni solid, light grey, dotted hatch
ax.add_patch(Rectangle((left, bot_band[0]), right - left,
                       bot_band[1] - bot_band[0], facecolor="0.88",
                       edgecolor="black", linewidth=0.9, hatch="...."))
ax.text((left + right) / 2, (bot_band[0] + bot_band[1]) / 2 + 0.04,
        "INNER CORE SOLID", ha="center", va="center", fontweight="bold",
        fontsize=8.5, color="black")
ax.text((left + right) / 2, (bot_band[0] + bot_band[1]) / 2 - 0.05,
        "Fe Ni, O 0.003, S 0.011 mol %", ha="center", va="center",
        fontsize=7.5, color="black")

# Freezing front line between the two reservoirs
front_y = 0.525
ax.plot([left, right], [front_y, front_y], color="black", lw=1.3, ls="-")
ax.text(right + 0.01, front_y, "freezing\nfront", ha="left", va="center",
        fontsize=6.8, color="black")

# Partition arrows downward through the front, with D_i labels
for xc, lab in [(0.45, "O  D 0.02"), (0.93, "S  D 0.75"),
                (1.41, "Ni  D 1.0")]:
    ar = FancyArrowPatch((xc, top_band[0] - 0.01), (xc, bot_band[1] + 0.01),
                         arrowstyle="-|>", mutation_scale=12, color="black",
                         lw=1.2, shrinkA=0, shrinkB=0)
    ax.add_patch(ar)
    ax.text(xc + 0.015, front_y, lab, ha="left", va="center", fontsize=6.6,
            color="black")

# One line outcome banner
ax.text(W / 2, 0.085,
        "Oxygen required in the outer core; nickel closes the inner core "
        "density", ha="center", va="center", fontsize=8.2, color="black",
        style="italic")

# Source note
fig.text(0.5, 0.018, SOURCE, ha="center", va="bottom", fontsize=7,
         color="black", wrap=True)

out = os.path.join(FIGDIR, "graphical_abstract.png")
fig.savefig(out, dpi=300)
plt.close(fig)
print("wrote", out, "ratio", round(W, 3))
