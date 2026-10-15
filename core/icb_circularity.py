"""Circularity audit of the oxygen requirement (referee response), on the production setup.

Test A: can an oxygen-free liquid satisfy the double constraint? Searched with the production
literature caps (C 1.0, Si 2.4, S 1.9 wt%) and, for comparison, with the caps lifted, which shows
what an oxygen-free core would have to carry. Test B: degrade the oxygen component off its Badro fit
toward iron; is oxygen still required within the caps?
Run: python icb_circularity.py
"""
import copy
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import full_space_rerun as F
import icb_partition_model as M

F._M = M
CAPS = {"C": 1.0, "Si": 2.4, "S": 1.9}
GRID = dict(O=np.round(np.arange(0, 0.205, 0.01), 3), Si=np.round(np.arange(0, 0.125, 0.01), 3),
            S=np.round(np.arange(0, 0.065, 0.005), 4), C=np.round(np.arange(0, 0.205, 0.01), 3))


def search(o_values, caps):
    F._CAPS_WT = dict(caps)
    M.base._phase_props.cache_clear()
    return F._sweep({**GRID, "O": np.asarray(o_values)}, M.T_ICB)


def describe(df):
    a = df[df.admissible]
    if a.empty:
        b = df.loc[df.rms_misfit_pct.idxmin()]
        w = F._wt_pct(b.liq_x_o, b.liq_x_si, b.liq_x_s, b.liq_x_c)
        fails = [k for k, t in F.TOL.items() if abs(b[k]) > t]
        return (f"none admissible; closest Si {w['Si']:.1f} S {w['S']:.1f} C {w['C']:.1f} wt%, "
                f"rms {b.rms_misfit_pct:.2f}%, fails {fails}")
    w = [F._wt_pct(r.liq_x_o, r.liq_x_si, r.liq_x_s, r.liq_x_c) for r in a.itertuples()]
    return (f"{len(a)} admissible; needs S {min(x['S'] for x in w):.1f}-{max(x['S'] for x in w):.1f}, "
            f"C {min(x['C'] for x in w):.1f}-{max(x['C'] for x in w):.1f}, "
            f"Si {min(x['Si'] for x in w):.1f}-{max(x['Si'] for x in w):.1f} wt%")


if __name__ == "__main__":
    print("=== TEST A: can an OXYGEN-FREE liquid satisfy the double constraint? ===")
    print("  within the literature caps:", describe(search([0.0], CAPS)))
    print("  caps lifted:               ", describe(search([0.0], {})))

    print("\n=== TEST B: degrade the oxygen calibration (within the caps); is oxygen still required? ===")
    liq = M.PP["liquid_outer_core"]
    orig = copy.deepcopy(liq)
    for scale in (1.0, 0.7, 0.5, 0.3, 0.0):
        for p in ("130", "329"):
            c, o = liq[f"components_at_{p}GPa"], orig[f"components_at_{p}GPa"]
            c["O"]["V"] = o["Fe"]["V"] - (o["Fe"]["V"] - o["O"]["V"]) * scale
            c["O"]["K"] = o["Fe"]["K"] - (o["Fe"]["K"] - o["O"]["K"]) * scale
        a = search(GRID["O"], CAPS)
        a = a[a.admissible]
        o_wt = [F._wt_pct(r.liq_x_o, r.liq_x_si, r.liq_x_s, r.liq_x_c)["O"] for r in a.itertuples()]
        req = len(a) > 0 and not (a.liq_x_o <= 0.03).any()
        rng = f"{min(o_wt):.1f}-{max(o_wt):.1f} wt%" if o_wt else "none"
        print(f"  O potency {scale:.1f}: admissible {len(a):3d}, O {rng}, oxygen required {req}")
    M.PP["liquid_outer_core"] = orig
