"""
Greyscale figure suite for the Fe Ni Si O C inner core partition model.

Four figures, strict greyscale house style (no color; series distinguished by
line style, greyscale shade, and diagonal hatching). All curves carry the
cosmochemical 5 at% nickel in the solid metal; a thin pure Fe (no nickel)
reference is shown in figures 1 and 2.

Run:  python3 fe_si_o_c_figures.py
"""
import os
import warnings

warnings.filterwarnings("ignore")  # suppress BurnMan autograd warnings

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import icb_partition_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures")
os.makedirs(FIGDIR, exist_ok=True)

T_K = 6000.0
X_NI = 0.05
X_S = 0.0
SOURCE = ("Source: BurnMan mgd3 aggregate with 5 at% nickel "
          "(icb_partition_model); PREM Dziewonski & Anderson 1981; "
          "Ni Vocadlo 2006, McDonough 2003. Greyscale per house figure standard.")

# Greyscale shades and line styles per series (no color anywhere)
PREM_C = "0.0"
SHADE = {"FeNi": "0.30", "C": "0.45", "Si": "0.60", "O": "0.20", "pureFe": "0.50"}
LS = {"FeNi": "-", "C": "--", "Si": "-.", "O": (0, (1, 1)), "pureFe": (0, (3, 1, 1, 1))}

P_GPA = np.linspace(330.0, 364.0, 35)


def _box(ax):
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_color("black")
        s.set_linewidth(0.8)
    ax.tick_params(direction="out", length=3, width=0.6, colors="black")
    ax.set_facecolor("white")
    ax.grid(False)


def _source(fig):
    fig.text(0.5, 0.012, SOURCE, ha="center", va="bottom", fontsize=7, color="black")


def _legend(ax, **kw):
    leg = ax.legend(frameon=True, fontsize=8, **kw)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    leg.get_frame().set_facecolor("white")
    return leg


# ------- aggregate profiles vs pressure -------
def _profile(xc, xsi, xo, x_ni):
    rho = np.empty_like(P_GPA)
    vphi = np.empty_like(P_GPA)
    for i, p in enumerate(P_GPA):
        r, _k, v = M.solid_props(xc, xsi, xo, X_S, p * 1e9, T_K, x_ni=x_ni)
        rho[i] = r
        vphi[i] = v
    return rho, vphi


PREM_RHO = np.array([M.base.prem_inner_core(p)[0] for p in P_GPA])
PREM_VPHI = np.array([M.base.prem_inner_core(p)[1] for p in P_GPA])

MW = {"Fe": 55.845, "Ni": 58.693, "O": 15.999, "Si": 28.0855, "C": 12.011}


def wt_pct(xc=0.0, xsi=0.0, xo=0.0, xni=X_NI):
    """Mole fractions (iron the balance) to wt%; the model is defined on mole fractions."""
    x = {"C": xc, "Si": xsi, "O": xo, "Ni": xni}
    x["Fe"] = 1.0 - sum(x.values())
    tot = sum(x[k] * MW[k] for k in x)
    return {k: 100 * x[k] * MW[k] / tot for k in x}


SERIES = [
    ("FeNi", f"Fe with {wt_pct()['Ni']:.1f} wt% Ni", 0.0, 0.0, 0.0, X_NI),
    ("C", f"Fe Ni + {wt_pct(xc=0.05)['C']:.1f} wt% C", 0.05, 0.0, 0.0, X_NI),
    ("Si", f"Fe Ni + {wt_pct(xsi=0.10)['Si']:.1f} wt% Si", 0.0, 0.10, 0.0, X_NI),
    ("O", f"Fe Ni + {wt_pct(xo=0.10)['O']:.1f} wt% O", 0.0, 0.0, 0.10, X_NI),
]
PROFILES = {key: _profile(xc, xsi, xo, xni)
            for key, _lbl, xc, xsi, xo, xni in SERIES}
PURE_RHO, PURE_VPHI = _profile(0.0, 0.0, 0.0, 0.0)


def fig1_density():
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(P_GPA, PREM_RHO, color=PREM_C, lw=2.6, ls="-", label="PREM", zorder=5)
    for key, lbl, *_ in SERIES:
        rho, _ = PROFILES[key]
        ax.plot(P_GPA, rho, color=SHADE[key], lw=1.6, ls=LS[key], label=lbl)
    ax.plot(P_GPA, PURE_RHO, color=SHADE["pureFe"], lw=0.9, ls=LS["pureFe"],
            label="pure Fe, no nickel")
    _box(ax)
    ax.set_xlabel("Pressure, GPa", fontsize=10)
    ax.set_ylabel("Density, g cm$^{-3}$", fontsize=10)
    _legend(ax, loc="upper left")
    ax.set_xlim(P_GPA.min(), P_GPA.max())
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(os.path.join(FIGDIR, "fig1_density_vs_P.png"), dpi=600, metadata={"Software": None})
    plt.close(fig)


def fig2_vphi():
    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.plot(P_GPA, PREM_VPHI, color=PREM_C, lw=2.6, ls="-", label="PREM", zorder=5)
    for key, lbl, *_ in SERIES:
        _, vphi = PROFILES[key]
        ax.plot(P_GPA, vphi, color=SHADE[key], lw=1.6, ls=LS[key], label=lbl)
    ax.plot(P_GPA, PURE_VPHI, color=SHADE["pureFe"], lw=0.9, ls=LS["pureFe"],
            label="pure Fe, no nickel")
    _box(ax)
    ax.set_xlabel("Pressure, GPa", fontsize=10)
    ax.set_ylabel("Bulk sound velocity, km s$^{-1}$", fontsize=10)
    _legend(ax, loc="upper left")
    ax.set_xlim(P_GPA.min(), P_GPA.max())
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(os.path.join(FIGDIR, "fig2_vphi_vs_P.png"), dpi=600, metadata={"Software": None})
    plt.close(fig)


# ------- misfit map -------
def _misfit_grid(xc_fixed, xsi_grid, xo_grid, xc_grid, mode):
    P = 330.0
    target = M.base.prem_inner_core(P)[0]
    if mode == "SiO":
        X, Y = np.meshgrid(xsi_grid, xo_grid)
        Z = np.full(X.shape, np.nan)
        for j in range(X.shape[0]):
            for i in range(X.shape[1]):
                r = M.solid_props(0.0, X[j, i], Y[j, i], 0.0, P * 1e9, T_K, x_ni=X_NI)
                if r is not None:
                    Z[j, i] = 100.0 * (r[0] - target) / target
        return X, Y, Z
    else:  # CO
        X, Y = np.meshgrid(xc_grid, xo_grid)
        Z = np.full(X.shape, np.nan)
        for j in range(X.shape[0]):
            for i in range(X.shape[1]):
                r = M.solid_props(X[j, i], 0.0, Y[j, i], 0.0, P * 1e9, T_K, x_ni=X_NI)
                if r is not None:
                    Z[j, i] = 100.0 * (r[0] - target) / target
        return X, Y, Z


def fig3_misfit_map():
    si = np.arange(0.0, 0.1001, 0.01)
    o = np.arange(0.0, 0.1001, 0.01)
    c = np.arange(0.0, 0.0501, 0.005)

    Xa, Ya, Za = _misfit_grid(0.0, si, o, c, "SiO")
    Xb, Yb, Zb = _misfit_grid(0.0, si, o, c, "CO")

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.0, 5.0))
    greys = plt.get_cmap("Greys")

    wa = [[wt_pct(xsi=x, xo=y) for x, y in zip(rx, ry)] for rx, ry in zip(Xa, Ya)]
    wb = [[wt_pct(xc=x, xo=y) for x, y in zip(rx, ry)] for rx, ry in zip(Xb, Yb)]
    XaW = np.array([[w["Si"] for w in r] for r in wa]); YaW = np.array([[w["O"] for w in r] for r in wa])
    XbW = np.array([[w["C"] for w in r] for r in wb]); YbW = np.array([[w["O"] for w in r] for r in wb])
    for ax, (X, Y, Z), xlab, header in (
        (axA, (XaW, YaW, Za), "Silicon, wt%", "(a)"),
        (axB, (XbW, YbW, Zb), "Carbon, wt%", "(b)"),
    ):
        zmin = np.nanmin(Z)
        zmax = np.nanmax(Z)
        levels = np.linspace(zmin, zmax, 13)
        cf = ax.contourf(X, Y, Z, levels=levels, cmap=greys)
        cl = ax.contour(X, Y, Z, levels=levels, colors="black", linewidths=0.6)
        ax.clabel(cl, inline=True, fontsize=7, fmt="%.1f")
        if zmin < 0.0 < zmax:
            zc = ax.contour(X, Y, Z, levels=[0.0], colors="black", linewidths=2.2)
            ax.clabel(zc, inline=True, fontsize=8, fmt="%.0f")
        ax.set_xlim(0.0, np.nanmin(X[:, -1]))   # mole-fraction grid mapped to wt%: crop to full coverage
        ax.set_ylim(0.0, np.nanmin(Y[-1, :]))
        _box(ax)
        ax.set_xlabel(xlab, fontsize=10)
        ax.set_ylabel("Oxygen, wt%", fontsize=10)
        ax.set_title(header, loc="left", fontsize=10, fontweight="bold")
        cb = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cb.outline.set_edgecolor("black")
        cb.outline.set_linewidth(0.8)
        cb.ax.tick_params(width=0.6, color="black", labelsize=8)
        cb.set_label("Density misfit vs PREM, %", fontsize=9)

    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    fig.savefig(os.path.join(FIGDIR, "fig3_misfit_map.png"), dpi=600, metadata={"Software": None})
    plt.close(fig)


# ------- carbon arm -------
def fig4_carbon_arm():
    xc = np.linspace(0.0, 0.05, 26)
    P = 330.0
    t_rho, _, t_vp, t_vs = M.base.prem_inner_core(P)
    t_vphi = M.base.prem_inner_core(P)[1]
    drho = np.empty_like(xc)
    dvphi = np.empty_like(xc)
    for i, c in enumerate(xc):
        r = M.solid_props(c, 0.0, 0.0, 0.0, P * 1e9, T_K, x_ni=X_NI)
        drho[i] = 100.0 * (r[0] - t_rho) / t_rho
        dvphi[i] = 100.0 * (r[2] - t_vphi) / t_vphi

    fig, ax = plt.subplots(figsize=(7.2, 5.0))
    ax.axhline(0.0, color="black", lw=0.9, ls=(0, (4, 3)), zorder=1)
    cw = np.array([wt_pct(xc=c)["C"] for c in xc])
    l1, = ax.plot(cw, drho, color="0.20", lw=1.8, ls="-",
                  label="Density misfit")
    l2, = ax.plot(cw, dvphi, color="0.50", lw=1.8, ls="--",
                  label="Bulk sound velocity misfit")
    # diagonal hatching fill between curve and zero to distinguish series
    ax.fill_between(cw, drho, 0.0, facecolor="none", edgecolor="0.20",
                    hatch="////", linewidth=0.0, alpha=0.6, zorder=0)
    _box(ax)
    ax.set_xlabel("Carbon, wt%", fontsize=10)
    ax.set_ylabel("Misfit vs PREM, %", fontsize=10)
    ax.set_xlim(0.0, cw[-1])
    _legend(ax, loc="best")
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    fig.savefig(os.path.join(FIGDIR, "fig4_carbon_arm.png"), dpi=600, metadata={"Software": None})
    plt.close(fig)


if __name__ == "__main__":
    fig1_density()
    fig2_vphi()
    fig3_misfit_map()
    fig4_carbon_arm()
    r = M.solid_props(0.0, 0.0, 0.0, 0.0, 329e9, T_K, x_ni=X_NI)[0]
    p = M.base.prem_inner_core(329)[0]
    print("pure Fe + Ni density misfit @329 GPa: %+.3f%%" % (100.0 * (r - p) / p))
    print("figures written to", FIGDIR)
