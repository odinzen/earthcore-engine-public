"""Two checks on the production setup that the paper quotes (output: production/checks.txt).

1. Carbon partition coefficient: production D_C 1.5 is a scenario value for carbon favouring the
   solid, so the reference-state search is repeated at D_C 1.0 and 2.0 and at Li 2019's 6.3e-5.
2. Oxygen component: the Badro-calibrated oxygen is replaced by the independent Xie 2025 component
   (liquid_refit_independent_data.py) and the search repeated.
All searches use the production caps (C 1.0, Si 2.4, S 1.9 wt%) and grid.
Run:  python production_checks.py
"""
import copy
import os
import warnings

import numpy as np

warnings.filterwarnings("ignore")
import full_space_rerun as F
import icb_partition_model as M
import liquid_refit_independent_data as R

F._M = M
F._CAPS_WT = {"C": 1.0, "Si": 2.4, "S": 1.9}
GRID = dict(O=np.round(np.arange(0, 0.205, 0.01), 3), Si=np.round(np.arange(0, 0.045, 0.01), 3),
            S=np.round(np.arange(0, 0.0305, 0.005), 4), C=np.round(np.arange(0, 0.045, 0.01), 3))
HERE = os.path.dirname(os.path.abspath(__file__))


def run(label):
    M.base._phase_props.cache_clear()
    df = F._sweep(GRID, M.T_ICB)
    a = df[df.admissible]
    if a.empty:
        return f"{label}: 0 admissible of {len(df)}"
    w = [F._wt_pct(r.liq_x_o, r.liq_x_si, r.liq_x_s, r.liq_x_c) for r in a.itertuples()]
    b = a.loc[a.rms_misfit_pct.idxmin()]
    wb = F._wt_pct(b.liq_x_o, b.liq_x_si, b.liq_x_s, b.liq_x_c)
    return (f"{label}: {len(a)} admissible of {len(df)}; O {min(x['O'] for x in w):.1f}-{max(x['O'] for x in w):.1f} wt%; "
            f"oxygen required {not (a.liq_x_o <= 0.03).any()}; best O {wb['O']:.1f} Si {wb['Si']:.1f} S {wb['S']:.1f} "
            f"C {wb['C']:.1f} wt%, rms {b.rms_misfit_pct:.2f}, jump {b.icb_jump_pct:.2f}")


if __name__ == "__main__":
    lines = []
    pc = M.PP["partition_coefficients"]["C"]
    d0 = pc["D"]
    for d in (d0, 1.0, 2.0, pc["D_scenario_li2019"]):
        pc["D"] = d
        lines.append(run(f"D_C {d:g}"))
    pc["D"] = d0
    orig = copy.deepcopy(M.PP["liquid_outer_core"])
    R.apply(("O",))
    lines.append(run("Xie 2025 oxygen component"))
    M.PP["liquid_outer_core"].update(orig)
    out = "\n".join(lines)
    print(out)
    open(os.path.join(HERE, "production", "checks.txt"), "w").write(out + "\n")
