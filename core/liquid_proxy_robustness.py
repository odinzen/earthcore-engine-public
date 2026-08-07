"""Liquid-proxy robustness of the oxygen-required conclusion.

The liquid outer-core model is a semi-empirical proxy: calibrated to Badro 2014 at
two pressures, interpolated linearly in P, temperature-independent, with the carbon
component set equal to oxygen and silicon an estimate. This checks whether the
headline result (oxygen required in the outer core) depends on those weak
assumptions, by perturbing them and re-running the Stage-2 admissibility search.
The search grid and tolerances mirror icb_partition_sweep.py; results are held in
memory so the committed sweep CSVs are never overwritten.

Run:  python liquid_proxy_robustness.py
"""
import copy
import warnings

warnings.filterwarnings("ignore")
import numpy as np

import icb_partition_model as M

X_O = np.round(np.arange(0.00, 0.205, 0.01), 3)
X_SI = np.round(np.arange(0.00, 0.105, 0.01), 3)
X_S = np.round(np.arange(0.00, 0.0605, 0.005), 4)
X_C = np.round(np.arange(0.00, 0.065, 0.02), 3)
TOL = dict(lr=1.0, lv=2.0, sr=1.5, sv=2.5)

_LIQ = M.PP["liquid_outer_core"]
_BASE = copy.deepcopy({"c130": _LIQ["components_at_130GPa"],
                       "c329": _LIQ["components_at_329GPa"]})


def reset():
    _LIQ["components_at_130GPa"] = copy.deepcopy(_BASE["c130"])
    _LIQ["components_at_329GPa"] = copy.deepcopy(_BASE["c329"])
    M._C1 = _LIQ["components_at_130GPa"]
    M._C2 = _LIQ["components_at_329GPa"]


def scale_component(el, fV=1.0, fK=1.0):
    for c in (_LIQ["components_at_130GPa"], _LIQ["components_at_329GPa"]):
        c[el]["V"] *= fV
        c[el]["K"] *= fK
    M._C1 = _LIQ["components_at_130GPa"]
    M._C2 = _LIQ["components_at_329GPa"]


def scale_all_density(fV=1.0):
    for el in ("Fe", "O", "S", "C", "Si", "Ni"):
        scale_component(el, fV=fV)


def search():
    pr_cmb = M.prem_outer(M.P_CMB)
    pr_icb = M.prem_outer(M.P_ICB)
    adm_o = []
    best = None
    best_rms = np.inf
    any_low_o = False
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
                    sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"],
                                       M.P_ICB * 1e9, M.T_ICB)
                    if sp is None:
                        continue
                    d = [100 * (lc[0] - pr_cmb[0]) / pr_cmb[0],
                         100 * (li[0] - pr_icb[0]) / pr_icb[0],
                         100 * (lc[1] - pr_cmb[1]) / pr_cmb[1],
                         100 * (li[1] - pr_icb[1]) / pr_icb[1],
                         100 * (sp[0] - M.PREM_IC["rho"]) / M.PREM_IC["rho"],
                         100 * (sp[2] - M.PREM_IC["vphi"]) / M.PREM_IC["vphi"]]
                    adm = (abs(d[0]) <= TOL["lr"] and abs(d[1]) <= TOL["lr"]
                           and abs(d[2]) <= TOL["lv"] and abs(d[3]) <= TOL["lv"]
                           and abs(d[4]) <= TOL["sr"] and abs(d[5]) <= TOL["sv"])
                    if not adm:
                        continue
                    adm_o.append(xo)
                    if xo <= 0.03:
                        any_low_o = True
                    rms = np.sqrt(np.mean(np.square(d)))
                    if rms < best_rms:
                        best_rms = rms
                        jump = 100 * (sp[0] - li[0]) / li[0]
                        best = (xo * 100, xs * 100, jump)
    if not adm_o:
        return dict(n=0, omin=None, omax=None, oxreq=True, best=None)
    return dict(n=len(adm_o), omin=min(adm_o) * 100, omax=max(adm_o) * 100,
                oxreq=(not any_low_o), best=best)


def show(tag, r):
    if r["n"] == 0:
        print(f"{tag:38s} NO admissible solution")
        return
    b = r["best"]
    print(f"{tag:38s} adm={r['n']:3d}  O={r['omin']:.0f}-{r['omax']:.0f} mol%  "
          f"O-required={r['oxreq']!s:5s}  best O={b[0]:.0f} S={b[1]:.1f} jump={b[2]:.1f}%")


def main():
    print("Liquid-proxy robustness: does 'oxygen required' survive perturbing the weak assumptions?\n")
    reset(); show("baseline (as published)", search())

    print("\n-- carbon component: currently assumed equal to oxygen --")
    reset(); scale_component("C", fV=1.30); show("  C molar volume +30% vs O", search())
    reset(); scale_component("C", fV=0.70); show("  C molar volume -30% vs O", search())

    print("\n-- silicon component: currently an interpolated estimate --")
    reset(); scale_component("Si", fV=1.30); show("  Si molar volume +30%", search())
    reset(); scale_component("Si", fV=0.70); show("  Si molar volume -30%", search())

    print("\n-- temperature-independence: shift whole-liquid density (thermal proxy) --")
    reset(); scale_all_density(fV=1.01); show("  liquid density -1% (hotter)", search())
    reset(); scale_all_density(fV=0.99); show("  liquid density +1% (cooler)", search())
    reset(); scale_all_density(fV=1.02); show("  liquid density -2% (hotter)", search())
    reset()


if __name__ == "__main__":
    main()
