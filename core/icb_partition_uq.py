"""Monte Carlo uncertainty quantification for the ICB partition double constraint.

Resamples partition coefficients, T_ICB, the pure Fe liquid CMB anchor, and the
solid EoS parameters; for each draw runs the admissibility search over outer core
O and S and records the best fit. Writes icb_partition_uq_results.csv and fig9_uq.png.
"""
import os
import warnings
import copy
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import icb_partition_model as M

HERE = os.path.dirname(os.path.abspath(__file__))
N_DRAWS = 400
RNG = np.random.default_rng(20260611)

SOLID_KEYS = ("hcp_fe", "feo", "fe7c3", "fe3s")
EL_KEYS = ("O", "S", "Si", "C", "Ni")

# Outer core search grid
X_O = np.round(np.arange(0.02, 0.241, 0.01), 4)
X_S = np.round(np.arange(0.00, 0.0401, 0.005), 4)

TOL_LIQ_RHO = 1.0
TOL_LIQ_VPHI = 2.0
TOL_SOL_RHO = 1.5
TOL_SOL_VPHI = 2.5

# Anchor central values for the pure Fe 130 GPa liquid
VPHI0 = 7.70
RHO0 = 10.70

SOURCE_NOTE = ("Source: Monte Carlo / sensitivity over EoS, partition, liquid anchor, "
               "and T_ICB (icb_partition_uq / icb_partition_sensitivity). "
               "Greyscale per house figure standard.")


def snapshot():
    """Deep copy the original mutable params we perturb each draw."""
    snap = {}
    for key in SOLID_KEYS:
        snap[key] = {"K_0": M.base._PHASES[key].params["K_0"],
                     "V_0": M.base._PHASES[key].params["V_0"]}
    liq_fe = M.PP["liquid_outer_core"]["components_at_130GPa"]["Fe"]
    snap["liq_fe"] = {"K": liq_fe["K"], "V": liq_fe["V"]}
    snap["D"] = {el: M.PP["partition_coefficients"][el]["D"] for el in EL_KEYS}
    return copy.deepcopy(snap)


def reset(snap):
    for key in SOLID_KEYS:
        M.base._PHASES[key].params["K_0"] = snap[key]["K_0"]
        M.base._PHASES[key].params["V_0"] = snap[key]["V_0"]
    liq_fe = M.PP["liquid_outer_core"]["components_at_130GPa"]["Fe"]
    liq_fe["K"] = snap["liq_fe"]["K"]
    liq_fe["V"] = snap["liq_fe"]["V"]
    for el in EL_KEYS:
        M.PP["partition_coefficients"][el]["D"] = snap["D"][el]
    M.base._phase_props.cache_clear()


def search(t_icb):
    """Admissibility search over outer core O and S at this draw's params.

    Returns the best fit row (min RMS misfit) and a flag for any admissible O<=0.03.
    """
    pr_cmb = M.prem_outer(M.P_CMB)
    pr_icb = M.prem_outer(M.P_ICB)
    best = None
    best_rms = np.inf
    any_adm = False
    any_low_o = False
    for xo in X_O:
        for xs in X_S:
            lc = M.liquid_props(0.0, xo, 0.0, xs, M.P_CMB)
            li = M.liquid_props(0.0, xo, 0.0, xs, M.P_ICB)
            if lc is None or li is None:
                continue
            xsol = M.partition_to_solid({"C": 0.0, "O": xo, "Si": 0.0, "S": xs})
            sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"],
                               M.P_ICB * 1e9, t_icb)
            if sp is None:
                continue
            d_lc_rho = 100 * (lc[0] - pr_cmb[0]) / pr_cmb[0]
            d_lc_v = 100 * (lc[1] - pr_cmb[1]) / pr_cmb[1]
            d_li_rho = 100 * (li[0] - pr_icb[0]) / pr_icb[0]
            d_li_v = 100 * (li[1] - pr_icb[1]) / pr_icb[1]
            d_s_rho = 100 * (sp[0] - M.PREM_IC["rho"]) / M.PREM_IC["rho"]
            d_s_v = 100 * (sp[2] - M.PREM_IC["vphi"]) / M.PREM_IC["vphi"]
            adm = (abs(d_lc_rho) <= TOL_LIQ_RHO and abs(d_li_rho) <= TOL_LIQ_RHO
                   and abs(d_lc_v) <= TOL_LIQ_VPHI and abs(d_li_v) <= TOL_LIQ_VPHI
                   and abs(d_s_rho) <= TOL_SOL_RHO and abs(d_s_v) <= TOL_SOL_VPHI)
            if not adm:
                continue
            any_adm = True
            if xo <= 0.03:
                any_low_o = True
            jump = 100 * (sp[0] - li[0]) / li[0]
            rms = np.sqrt(np.mean(np.square(
                [d_lc_rho, d_li_rho, d_lc_v, d_li_v, d_s_rho, d_s_v])))
            if rms < best_rms:
                best_rms = rms
                best = dict(oc_o=xo, oc_s=xs, sol_o=xsol["O"], sol_s=xsol["S"],
                            d_liq_rho_cmb=d_lc_rho, d_liq_vphi_cmb=d_lc_v,
                            d_liq_rho_icb=d_li_rho, d_liq_vphi_icb=d_li_v,
                            d_sol_rho_icb=d_s_rho, d_sol_vphi_icb=d_s_v,
                            icb_jump_pct=jump, rms_misfit_pct=rms)
    return best, any_adm, any_low_o


def atom_o_to_wt(x_o):
    """Convert mol fraction O (rest Fe) to wt% O."""
    mo = M.M_ATOM["O"]
    mfe = M.M_ATOM["Fe"]
    return 100.0 * x_o * mo / (x_o * mo + (1.0 - x_o) * mfe)


def atom_s_to_wt(x_s):
    ms = M.M_ATOM["S"]
    mfe = M.M_ATOM["Fe"]
    return 100.0 * x_s * ms / (x_s * ms + (1.0 - x_s) * mfe)


def main():
    snap = snapshot()
    rows = []
    for i in range(N_DRAWS):
        reset(snap)

        # Resample partition coefficients
        d_o = 0.02 * 2.0 ** RNG.uniform(-1.0, 1.0)
        d_s = max(1e-6, RNG.normal(0.75, 0.07))
        d_si = RNG.uniform(0.67, 1.0)
        d_ni = max(1e-6, RNG.normal(1.0, 0.1))
        d_c = 6.3e-5

        M.PP["partition_coefficients"]["O"]["D"] = d_o
        M.PP["partition_coefficients"]["S"]["D"] = d_s
        M.PP["partition_coefficients"]["Si"]["D"] = d_si
        M.PP["partition_coefficients"]["Ni"]["D"] = d_ni
        M.PP["partition_coefficients"]["C"]["D"] = d_c

        t_icb = RNG.uniform(5500.0, 6200.0)

        # Liquid pure Fe CMB anchor perturbation -> scale 130 GPa Fe component
        vphi = RNG.normal(VPHI0, 0.20)
        rho = RNG.normal(RHO0, 0.20)
        liq_fe = M.PP["liquid_outer_core"]["components_at_130GPa"]["Fe"]
        liq_fe["K"] = snap["liq_fe"]["K"] * (vphi / VPHI0) ** 2 * (rho / RHO0)
        liq_fe["V"] = snap["liq_fe"]["V"] * (RHO0 / rho)

        # Solid EoS perturbation
        for key in SOLID_KEYS:
            fk = RNG.normal(1.0, 0.02)
            fv = RNG.normal(1.0, 0.02)
            M.base._PHASES[key].params["K_0"] = snap[key]["K_0"] * fk
            M.base._PHASES[key].params["V_0"] = snap[key]["V_0"] * fv
        # Fe3S K_0 additional absolute +/- 7 GPa
        M.base._PHASES["fe3s"].params["K_0"] += RNG.normal(0.0, 7.0e9)
        M.base._phase_props.cache_clear()

        best, any_adm, any_low_o = search(t_icb)

        row = dict(draw=i, d_o=d_o, d_s=d_s, d_si=d_si, d_ni=d_ni, d_c=d_c,
                   t_icb=t_icb, liq_vphi=vphi, liq_rho=rho,
                   admissible=bool(any_adm),
                   oxygen_not_required=bool(any_low_o))
        if best is not None:
            row.update(best)
            row["best_oc_o_wt"] = atom_o_to_wt(best["oc_o"])
            row["best_oc_s_wt"] = atom_s_to_wt(best["oc_s"])
            row["best_sol_s_wt"] = atom_s_to_wt(best["sol_s"])
        else:
            for k in ("oc_o", "oc_s", "sol_o", "sol_s", "d_liq_rho_cmb",
                      "d_liq_vphi_cmb", "d_liq_rho_icb", "d_liq_vphi_icb",
                      "d_sol_rho_icb", "d_sol_vphi_icb", "icb_jump_pct",
                      "rms_misfit_pct", "best_oc_o_wt", "best_oc_s_wt",
                      "best_sol_s_wt"):
                row[k] = np.nan
        rows.append(row)
        reset(snap)

    df = pd.DataFrame(rows)
    out_csv = os.path.join(HERE, "icb_partition_uq_results.csv")
    df.to_csv(out_csv, index=False)

    adm = df[df.admissible]
    frac_adm = len(adm) / len(df)
    # Oxygen required = admissible draw with NO admissible O<=3 mol%
    frac_o_required = (adm.oxygen_not_required == False).mean() if len(adm) else float("nan")

    def pctl(s):
        s = s.dropna()
        return (np.median(s), np.percentile(s, 16), np.percentile(s, 84))

    o_mol = pctl(adm.oc_o * 100.0)
    o_wt = pctl(adm.best_oc_o_wt)
    s_oc_mol = pctl(adm.oc_s * 100.0)
    s_ic_mol = pctl(adm.sol_s * 100.0)
    jump = pctl(adm.icb_jump_pct)

    print("=" * 64)
    print("MONTE CARLO UNCERTAINTY QUANTIFICATION  (N=%d draws)" % N_DRAWS)
    print("=" * 64)
    print(f"Fraction of draws with >=1 admissible solution : {frac_adm:.3f}")
    print(f"Fraction of admissible draws requiring oxygen   : {frac_o_required:.3f}")
    print("-" * 64)
    print("Quantity                 median   [16th , 84th ]")
    print(f"Outer core O (mol%)      {o_mol[0]:6.1f}   [{o_mol[1]:.1f} , {o_mol[2]:.1f}]")
    print(f"Outer core O (wt%)       {o_wt[0]:6.1f}   [{o_wt[1]:.1f} , {o_wt[2]:.1f}]")
    print(f"Outer core S (mol%)      {s_oc_mol[0]:6.1f}   [{s_oc_mol[1]:.1f} , {s_oc_mol[2]:.1f}]")
    print(f"Inner core S (mol%)      {s_ic_mol[0]:6.1f}   [{s_ic_mol[1]:.1f} , {s_ic_mol[2]:.1f}]")
    print(f"ICB density jump (%)     {jump[0]:6.1f}   [{jump[1]:.1f} , {jump[2]:.1f}]")
    print("=" * 64)

    make_fig(adm, o_mol, jump, out_csv)
    print(f"Wrote {out_csv}")


def style_axes(ax):
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color("black")
        sp.set_linewidth(0.8)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.tick_params(colors="black", labelsize=8)


def make_fig(adm, o_mol, jump, out_csv):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10, 4.2))
    fig.patch.set_facecolor("white")

    # PANEL A: best fit outer core O (mol%)
    o_vals = (adm.oc_o * 100.0).dropna()
    axA.hist(o_vals, bins=np.arange(1.5, 25.5, 1.0), color="0.65",
             edgecolor="black", linewidth=0.6)
    axA.axvline(o_mol[0], color="black", linewidth=1.6, label="median")
    axA.axvline(o_mol[1], color="0.35", linewidth=1.0, linestyle="--",
                label="16th to 84th")
    axA.axvline(o_mol[2], color="0.35", linewidth=1.0, linestyle="--")
    axA.set_title("PANEL A", fontsize=10, fontweight="bold")
    axA.set_xlabel("BEST FIT OUTER CORE O, mol percent")
    axA.set_ylabel("NUMBER OF DRAWS, count")
    leg = axA.legend(fontsize=8, frameon=True)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    style_axes(axA)

    # PANEL B: ICB density jump (%) with observed 4.5 to 4.8 band
    j_vals = adm.icb_jump_pct.dropna()
    axB.hist(j_vals, bins=np.arange(3.5, 7.6, 0.25), color="0.65",
             edgecolor="black", linewidth=0.6)
    axB.axvspan(4.5, 4.8, color="0.85", alpha=1.0, zorder=0,
                label="observed 4.5 to 4.8")
    axB.axvline(jump[0], color="black", linewidth=1.6, label="median")
    axB.set_title("PANEL B", fontsize=10, fontweight="bold")
    axB.set_xlabel("ICB DENSITY JUMP, percent")
    axB.set_ylabel("NUMBER OF DRAWS, count")
    leg = axB.legend(fontsize=8, frameon=True)
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    style_axes(axB)

    fig.suptitle("MONTE CARLO UNCERTAINTY OF THE ICB PARTITION DOUBLE CONSTRAINT",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.005, SOURCE_NOTE, ha="center", va="bottom", fontsize=7,
             color="black")
    fig.tight_layout(rect=[0, 0.045, 1, 0.95])
    out_png = os.path.join(os.path.dirname(out_csv), "fig9_uq.png")
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
