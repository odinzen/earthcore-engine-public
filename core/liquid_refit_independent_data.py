"""Refit of the liquid outer-core components to INDEPENDENT verified data, and the
resulting robustness of the oxygen-required conclusion.

Motivation: the production proxy calibrates all liquid components to Badro 2014 alone,
at two pressures, with carbon assumed equal to oxygen and silicon an estimate. This
re-anchors the silicon and oxygen components to independent first-principles datasets
and re-runs the Stage-2 sweep, to test whether "oxygen required, silicon excluded"
survives independent liquid data (both volume AND modulus, since silicon is
discriminated by velocity).

Sources (verified vs Crossref; component V,K derived here, not read from a table):
  Fe-Si  Jin, Xie, Fu & Belonoshko 2026, Minerals 16, 576 (10.3390/min16060576).
         Supplementary "FP-MD data.xls" P-V-T grid; component V_Si,K_Si isolated by a
         BM3 fit per isotherm + ideal-mixing decomposition over 0-10.3 wt% Si.
  Fe-O   Xie, Fu & Belonoshko 2025, Geosci. Front. 16, 101847 (10.1016/j.gsf.2024.101847).
         Supplementary mmc2.xlsx P-V-T grid; V_O,K_O by the same method over 0-9.9 wt% O.
  Fe-C   Nakajima et al. 2015, Nat. Commun. 6, 8942 (10.1038/ncomms9942) for the C
         adiabatic modulus (K_S0=110 GPa, K'_S=5.14).
Method note: FP-MD (DFT) pressures carry a ~+17 GPa correction to the experiment scale
(Jin Table S1); K_T from -V dP/dV, adiabatic K_S = K_T*(1+aGT) ~ 1.10 K_T. That factor
is validated in-place: the derived Fe component reproduces the model's Fe (V,K) to <2%
at both anchors, so the decomposition and pressure scale are consistent.

Derived component (V cm3/mol-atom, K_S GPa) at the 130 / 329 GPa anchors:
  Si : (5.709, 264) / (4.301, 1688)   [model uses (4.250, 623)/(4.000, 1250)]
  O  : (3.459, 505) / (2.589, 1137)   [model uses (3.289, 613)/(2.487, 1108); ~5-10% of Badro]

Result: with these components, oxygen is required and the silicon low-oxygen escape
route is velocity-excluded, robust up to +10% on K_Si (opens at +15%). The independent
oxygen (Xie) validates the Badro oxygen; the ~6% density jump persists (not a proxy
artifact). See journal 2026-07-31.

Run:  python liquid_refit_independent_data.py
"""
import copy
import warnings
warnings.filterwarnings("ignore")
import numpy as np

import icb_partition_model as M

DATA = {
    "Si": {130: (5.709, 264.0), 329: (4.301, 1688.0)},
    "O":  {130: (3.459, 505.0), 329: (2.589, 1137.0)},
}
X_O = np.round(np.arange(0.00, 0.205, 0.01), 3)
X_SI = np.round(np.arange(0.00, 0.105, 0.01), 3)
X_S = np.round(np.arange(0.00, 0.0605, 0.005), 4)
X_C = np.round(np.arange(0.00, 0.065, 0.02), 3)
TOL = dict(lr=1.0, lv=2.0, sr=1.5, sv=2.5)
_LIQ = M.PP["liquid_outer_core"]
_BASE = copy.deepcopy({130: _LIQ["components_at_130GPa"], 329: _LIQ["components_at_329GPa"]})


def reset():
    _LIQ["components_at_130GPa"].update(copy.deepcopy(_BASE[130]))
    _LIQ["components_at_329GPa"].update(copy.deepcopy(_BASE[329]))


def apply(elements, f_ksi=1.0):
    for el in elements:
        for P in (130, 329):
            v, k = DATA[el][P]
            if el == "Si":
                k *= f_ksi
            _LIQ[f"components_at_{P}GPa"][el]["V"] = v
            _LIQ[f"components_at_{P}GPa"][el]["K"] = k


def search():
    pr_cmb = M.prem_outer(M.P_CMB); pr_icb = M.prem_outer(M.P_ICB)
    adm_o = []; best = None; best_rms = np.inf; low_o = 0
    for xo in X_O:
        for xsi in X_SI:
            for xs in X_S:
                for xc in X_C:
                    if xo + xsi + xs + xc >= 0.6:
                        continue
                    lc = M.liquid_props(xc, xo, xsi, xs, M.P_CMB)
                    li = M.liquid_props(xc, xo, xsi, xs, M.P_ICB)
                    if lc is None or li is None:
                        continue
                    xsol = M.partition_to_solid({"C": xc, "O": xo, "Si": xsi, "S": xs})
                    sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"], M.P_ICB * 1e9, M.T_ICB)
                    if sp is None:
                        continue
                    d = [100 * (lc[0] - pr_cmb[0]) / pr_cmb[0], 100 * (li[0] - pr_icb[0]) / pr_icb[0],
                         100 * (lc[1] - pr_cmb[1]) / pr_cmb[1], 100 * (li[1] - pr_icb[1]) / pr_icb[1],
                         100 * (sp[0] - M.PREM_IC["rho"]) / M.PREM_IC["rho"], 100 * (sp[2] - M.PREM_IC["vphi"]) / M.PREM_IC["vphi"]]
                    if not (abs(d[0]) <= TOL["lr"] and abs(d[1]) <= TOL["lr"] and abs(d[2]) <= TOL["lv"]
                            and abs(d[3]) <= TOL["lv"] and abs(d[4]) <= TOL["sr"] and abs(d[5]) <= TOL["sv"]):
                        continue
                    adm_o.append(xo)
                    if xo <= 0.03:
                        low_o += 1
                    rms = np.sqrt(np.mean(np.square(d)))
                    if rms < best_rms:
                        best_rms = rms
                        best = (xo * 100, xs * 100, 100 * (sp[0] - li[0]) / li[0])
    if not adm_o:
        return "NO admissible solution"
    return (f"adm={len(adm_o):3d}  O={min(adm_o)*100:.0f}-{max(adm_o)*100:.0f} mol%  "
            f"O-required={(low_o == 0)!s:5s}  best O={best[0]:.0f} S={best[1]:.1f} jump={best[2]:.1f}%")


def main():
    print("Liquid refit against independent verified data (V and adiabatic K):\n")
    reset(); print(f"  baseline (Badro proxy)          {search()}")
    reset(); apply(["O"]); print(f"  O  -> Xie 2025                  {search()}")
    reset(); apply(["Si"]); print(f"  Si -> Jin 2026                  {search()}")
    reset(); apply(["Si", "O"]); print(f"  Si+O -> independent data        {search()}")
    print("\n  silicon-modulus robustness of the exclusion (K_Si scaled):")
    for f in (0.90, 1.00, 1.10, 1.15, 1.20):
        reset(); apply(["Si", "O"], f_ksi=f)
        print(f"    K_Si x{f:.2f}: {search()}")
    reset()


if __name__ == "__main__":
    main()
