"""Four-element rerun of the Stage-2 search, sensitivity, and Monte Carlo.

The production icb_partition_uq.py and icb_partition_sensitivity.py search outer core
O and S only (Si = C = 0), and the main sweep stops carbon at 6 mol%. Here every
search covers O, Si, S and C, with carbon extended to 24 mol%, so the oxygen
requirement is tested against every light element the model carries. The 400 Monte
Carlo draws reproduce icb_partition_uq.py exactly (same seed, same call order; checked
against icb_partition_uq_results.csv before anything runs).

Outputs go to rerun_full_space/ and nothing in the production CSVs or figures changes.

Run:  python full_space_rerun.py [--workers N] [--smoke]
"""
import argparse
import json
import os
import sys
import time
import warnings
from multiprocessing import Pool

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "rerun_full_space")

TOL = dict(lc_rho=1.0, li_rho=1.0, lc_v=2.0, li_v=2.0, s_rho=1.5, s_v=2.5)
LOW_O = 0.03  # same "oxygen not required" cut as icb_partition_uq.py

# Main sweep and sensitivity: the production steps, carbon and oxygen extended.
GRID_MAIN = dict(O=np.round(np.arange(0.00, 0.245, 0.01), 3),
                 Si=np.round(np.arange(0.00, 0.105, 0.01), 3),
                 S=np.round(np.arange(0.00, 0.0605, 0.005), 4),
                 C=np.round(np.arange(0.00, 0.245, 0.02), 3))
# Monte Carlo: the production O and S steps, Si and C added at 2 mol%.
GRID_UQ = dict(O=np.round(np.arange(0.00, 0.245, 0.01), 3),
               Si=np.round(np.arange(0.00, 0.105, 0.02), 3),
               S=np.round(np.arange(0.00, 0.0401, 0.005), 4),
               C=np.round(np.arange(0.00, 0.245, 0.02), 3))

SOLID_KEYS = ("hcp_fe", "feo", "fe7c3", "fe3s")
EL_KEYS = ("O", "S", "Si", "C", "Ni")
VPHI0, RHO0 = 7.70, 10.70

_M = None
_CAPS_WT = {}   # e.g. {"C": 1.5, "Si": 2.4, "S": 1.9}: compositions above any cap are not evaluated
MW = {"Fe": 55.845, "Ni": 58.693, "O": 15.999, "Si": 28.0855, "S": 32.06, "C": 12.011}


def _wt_pct(xo, xsi, xs, xc, x_ni=0.05):
    x = {"O": xo, "Si": xsi, "S": xs, "C": xc, "Ni": x_ni}
    x["Fe"] = 1.0 - sum(x.values())
    tot = sum(x[k] * MW[k] for k in x)
    return {k: 100 * x[k] * MW[k] / tot for k in ("O", "Si", "S", "C")}


def _init(feo_phase=None, carbon="equal_to_O", independent=False, caps=None):
    global _M, _CAPS_WT
    _CAPS_WT = dict(caps or {})
    if feo_phase:
        os.environ["FEO_PHASE"] = feo_phase
    sys.path.insert(0, HERE)
    import icb_partition_model as M
    if carbon == "badro":  # liquid_carbon_badro.json, from Badro 2014 Fig. 1
        comp = json.load(open(os.path.join(HERE, "liquid_carbon_badro.json")))["components"]
        for p in ("130", "329"):
            M.PP["liquid_outer_core"][f"components_at_{p}GPa"]["C"].update(comp[f"{p}GPa"])
    if independent:  # Si (Jin 2026) and O (Xie 2025) components of liquid_refit_independent_data.py
        import liquid_refit_independent_data as R
        R.apply(("Si", "O"))
    _M = M


def _eval(xc, xo, xsi, xs, t_icb, pr_cmb, pr_icb):
    M = _M
    lc = M.liquid_props(xc, xo, xsi, xs, M.P_CMB)
    li = M.liquid_props(xc, xo, xsi, xs, M.P_ICB)
    if lc is None or li is None:
        return None
    xsol = M.partition_to_solid({"C": xc, "O": xo, "Si": xsi, "S": xs})
    sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"], M.P_ICB * 1e9, t_icb)
    if sp is None:
        return None
    d = dict(lc_rho=100 * (lc[0] - pr_cmb[0]) / pr_cmb[0], lc_v=100 * (lc[1] - pr_cmb[1]) / pr_cmb[1],
             li_rho=100 * (li[0] - pr_icb[0]) / pr_icb[0], li_v=100 * (li[1] - pr_icb[1]) / pr_icb[1],
             s_rho=100 * (sp[0] - M.PREM_IC["rho"]) / M.PREM_IC["rho"],
             s_v=100 * (sp[2] - M.PREM_IC["vphi"]) / M.PREM_IC["vphi"])
    d["admissible"] = all(abs(d[k]) <= TOL[k] for k in TOL)
    d["rms_misfit_pct"] = float(np.sqrt(np.mean([d[k] ** 2 for k in TOL])))
    d["icb_jump_pct"] = 100 * (sp[0] - li[0]) / li[0]
    d.update(liq_x_c=xc, liq_x_o=xo, liq_x_si=xsi, liq_x_s=xs,
             sol_x_c=xsol["C"], sol_x_o=xsol["O"], sol_x_si=xsol["Si"], sol_x_s=xsol["S"])
    return d


def _sweep(grid, t_icb):
    M = _M
    pr_cmb, pr_icb = M.prem_outer(M.P_CMB), M.prem_outer(M.P_ICB)
    rows = []
    for xo in grid["O"]:
        for xsi in grid["Si"]:
            for xs in grid["S"]:
                for xc in grid["C"]:
                    if xo + xsi + xs + xc >= 0.5:
                        continue
                    if _CAPS_WT:
                        w = _wt_pct(xo, xsi, xs, xc)
                        if any(w[el] > cap + 1e-9 for el, cap in _CAPS_WT.items()):
                            continue
                    d = _eval(xc, xo, xsi, xs, t_icb, pr_cmb, pr_icb)
                    if d is not None:
                        rows.append(d)
    return pd.DataFrame(rows)


def summarize(df):
    """What the paper quotes from a search: the admissible box and the oxygen test."""
    a = df[df.admissible]
    if a.empty:
        return dict(n_admissible=0, admissible=False)
    best = a.loc[a.rms_misfit_pct.idxmin()]
    low = a[a.liq_x_o <= LOW_O]
    out = dict(n_tested=len(df), n_admissible=len(a), admissible=True,
               adm_o_min=a.liq_x_o.min(), adm_o_max=a.liq_x_o.max(),
               adm_si_max=a.liq_x_si.max(), adm_c_max=a.liq_x_c.max(),
               adm_s_min=a.liq_x_s.min(), adm_s_max=a.liq_x_s.max(),
               n_adm_with_si=int((a.liq_x_si > 0).sum()), n_adm_with_c=int((a.liq_x_c > 0).sum()),
               n_adm_low_o=len(low), oxygen_required=bool(low.empty),
               adm_o_plus_c_min=(a.liq_x_o + a.liq_x_c).min(),
               adm_o_plus_c_max=(a.liq_x_o + a.liq_x_c).max(),
               best_o=best.liq_x_o, best_si=best.liq_x_si, best_s=best.liq_x_s, best_c=best.liq_x_c,
               best_rms=best.rms_misfit_pct, best_jump=best.icb_jump_pct,
               best_sol_o=best.sol_x_o, best_sol_s=best.sol_x_s,
               jump_min=a.icb_jump_pct.min(), jump_max=a.icb_jump_pct.max())
    if not low.empty:
        b = low.loc[low.rms_misfit_pct.idxmin()]
        out.update(best_low_o=b.liq_x_o, best_low_o_c=b.liq_x_c, best_low_o_si=b.liq_x_si,
                   best_low_o_s=b.liq_x_s, best_low_o_rms=b.rms_misfit_pct)
    return out


# ---------------------------------------------------------------- stage tasks
def task_sweep(args):
    label, t_icb, d_si, grid = args
    M = _M
    if d_si is not None:
        M.PP["partition_coefficients"]["Si"]["D"] = d_si
    M.base._phase_props.cache_clear()
    df = _sweep(grid, t_icb)
    return label, df


def draw_params(n):
    """Replay the RNG call order of icb_partition_uq.main() (seed 20260611)."""
    rng = np.random.default_rng(20260611)
    draws = []
    for i in range(n):
        p = dict(draw=i)
        p["d_o"] = 0.02 * 2.0 ** rng.uniform(-1.0, 1.0)
        p["d_s"] = max(1e-6, rng.normal(0.75, 0.07))
        p["d_si"] = rng.uniform(0.67, 1.0)
        p["d_ni"] = max(1e-6, rng.normal(1.0, 0.1))
        p["d_c"] = 6.3e-5
        p["t_icb"] = rng.uniform(5500.0, 6200.0)
        p["liq_vphi"] = rng.normal(VPHI0, 0.20)
        p["liq_rho"] = rng.normal(RHO0, 0.20)
        p["solid"] = {k: (rng.normal(1.0, 0.02), rng.normal(1.0, 0.02)) for k in SOLID_KEYS}
        p["fe3s_dk"] = rng.normal(0.0, 7.0e9)
        draws.append(p)
    return draws


def task_draw(item):
    p, grid = item
    M = _M
    base = getattr(task_draw, "_snap", None)
    if base is None:
        base = dict(solid={k: (M.base._PHASES[k].params["K_0"], M.base._PHASES[k].params["V_0"])
                           for k in SOLID_KEYS},
                    liq_fe=dict(M.PP["liquid_outer_core"]["components_at_130GPa"]["Fe"]),
                    D={el: M.PP["partition_coefficients"][el]["D"] for el in EL_KEYS})
        task_draw._snap = base
    for el, key in (("O", "d_o"), ("S", "d_s"), ("Si", "d_si"), ("Ni", "d_ni"), ("C", "d_c")):
        M.PP["partition_coefficients"][el]["D"] = p[key]
    liq_fe = M.PP["liquid_outer_core"]["components_at_130GPa"]["Fe"]
    liq_fe["K"] = base["liq_fe"]["K"] * (p["liq_vphi"] / VPHI0) ** 2 * (p["liq_rho"] / RHO0)
    liq_fe["V"] = base["liq_fe"]["V"] * (RHO0 / p["liq_rho"])
    for k in SOLID_KEYS:
        fk, fv = p["solid"][k]
        M.base._PHASES[k].params["K_0"] = base["solid"][k][0] * fk
        M.base._PHASES[k].params["V_0"] = base["solid"][k][1] * fv
    M.base._PHASES["fe3s"].params["K_0"] += p["fe3s_dk"]
    M.base._phase_props.cache_clear()

    df = _sweep(grid, p["t_icb"])
    s = summarize(df)
    # The production Monte Carlo searched O and S only; its answer must reappear here.
    os_only = summarize(df[(df.liq_x_si == 0) & (df.liq_x_c == 0) & (df.liq_x_o >= 0.02)])
    s["os_only_best_o"] = os_only.get("best_o", np.nan)
    s["os_only_best_s"] = os_only.get("best_s", np.nan)

    for k in SOLID_KEYS:  # restore so the next draw in this worker starts clean
        M.base._PHASES[k].params["K_0"], M.base._PHASES[k].params["V_0"] = base["solid"][k]
    liq_fe.update(base["liq_fe"])
    for el in EL_KEYS:
        M.PP["partition_coefficients"][el]["D"] = base["D"][el]
    M.base._phase_props.cache_clear()

    row = {k: v for k, v in p.items() if k not in ("solid", "fe3s_dk")}
    row.update(s)
    return row


def write_best_fit(df, caps, out):
    """The reference-state best fit that the figures, velocities and text quote."""
    a = df[df.admissible]
    r = a.loc[a.rms_misfit_pct.idxmin()]
    liq = {el: float(r[f"liq_x_{el.lower()}"]) for el in ("O", "Si", "S", "C")}
    sol = {el: float(r[f"sol_x_{el.lower()}"]) for el in ("O", "Si", "S", "C")}
    best = {"_README": "Minimum-rms admissible liquid of production/sweep_baseline_5800K.csv (written by "
                       "full_space_rerun.py --production). Mole fractions; nickel 0.05 in both reservoirs.",
            "caps_wt_pct": caps, "liquid_x": liq, "solid_x": sol,
            "liquid_wt_pct": _wt_pct(liq["O"], liq["Si"], liq["S"], liq["C"]),
            "solid_wt_pct": _wt_pct(sol["O"], sol["Si"], sol["S"], sol["C"]),
            "misfit_pct": {k: float(r[k]) for k in TOL}, "rms_misfit_pct": float(r.rms_misfit_pct),
            "icb_jump_pct": float(r.icb_jump_pct), "n_admissible": int(len(a)), "n_tested": int(len(df))}
    json.dump(best, open(os.path.join(out, "best_fit.json"), "w"), indent=1)


def check_replay(draws):
    ref = pd.read_csv(os.path.join(HERE, "icb_partition_uq_results.csv"))
    got = pd.DataFrame([{k: d[k] for k in ("d_o", "d_s", "d_si", "d_ni", "t_icb", "liq_vphi", "liq_rho")}
                        for d in draws])
    n = min(len(ref), len(got))
    diff = np.abs(got.iloc[:n].values - ref[got.columns].iloc[:n].values).max()
    if diff > 1e-9:
        raise SystemExit(f"draw replay does not match icb_partition_uq_results.csv (max diff {diff:g})")
    print(f"draw replay matches the production Monte Carlo ({n} draws, max diff {diff:.1e})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--smoke", action="store_true", help="2 draws, coarse grids, separate output dir")
    ap.add_argument("--carbon", choices=("equal_to_O", "badro"), default="equal_to_O",
                    help="liquid carbon component: the production assumption or Badro 2014 Fig. 1")
    ap.add_argument("--independent", action="store_true",
                    help="use the independent Si (Jin 2026) and O (Xie 2025) liquid components")
    ap.add_argument("--anchor-sd", type=float, default=None,
                    help="relative 1-sigma of the CMB pure liquid Fe anchor (rho and vphi), e.g. 0.01 for "
                         "Kuwayama 2020; the same random draws are rescaled, so the replay still holds")
    ap.add_argument("--caps", default=None,
                    help='literature caps in wt%%, e.g. "C=1.5,Si=2.4,S=1.9"; grids shrink to cover them')
    ap.add_argument("--production", action="store_true",
                    help="the paper's setup: Kuwayama 1%% anchor, caps C=1.0,Si=2.4,S=1.9 wt%%, output to production/, "
                         "and best_fit.json for the figure and velocity scripts")
    ap.add_argument("--legacy-checks", action="store_true",
                    help="regression checks against the pre-2026-10-14 production CSVs (old model only)")
    ap.add_argument("--c-max", type=float, default=None,
                    help="cap liquid carbon (mol fraction), e.g. 0.04 for a literature carbon bound")
    args = ap.parse_args()

    global OUT, GRID_MAIN, GRID_UQ
    n_draws = 400
    if args.production:
        args.anchor_sd, args.caps = 0.01, "C=1.0,Si=2.4,S=1.9"
    if args.smoke:
        OUT = os.path.join(HERE, "rerun_full_space_smoke")
        GRID_MAIN = {k: v[::3] for k, v in GRID_MAIN.items()}
        GRID_UQ = {k: v[::3] for k, v in GRID_UQ.items()}
        n_draws = 2
    if args.carbon == "badro":
        OUT += "_carbon_badro"
    if args.independent:
        OUT += "_indep"
    if args.anchor_sd is not None:
        OUT += f"_anchor{round(args.anchor_sd * 1000):d}permil"
    if args.c_max is not None:
        OUT += f"_cmax{round(args.c_max * 100):d}"
        GRID_MAIN = {**GRID_MAIN, "C": GRID_MAIN["C"][GRID_MAIN["C"] <= args.c_max + 1e-9]}
        GRID_UQ = {**GRID_UQ, "C": GRID_UQ["C"][GRID_UQ["C"] <= args.c_max + 1e-9]}
    caps = {}
    if args.caps:
        caps = {k.strip(): float(v) for k, v in (kv.split("=") for kv in args.caps.split(","))}
        OUT += "_caps_" + "_".join(f"{k}{v:g}" for k, v in caps.items())
        # mol-fraction ceilings comfortably above each wt% cap, then the exact cap filters points
        top = {"C": caps.get("C", 24) * 4.6 / 100, "Si": caps.get("Si", 10) * 1.95 / 100, "S": caps.get("S", 6) * 1.75 / 100}
        GRID_MAIN = dict(O=np.round(np.arange(0, 0.205, 0.01), 3),
                         Si=np.round(np.arange(0, min(top["Si"], 0.10) + 1e-9, 0.01), 3),
                         S=np.round(np.arange(0, min(top["S"], 0.06) + 1e-9, 0.005), 4),
                         C=np.round(np.arange(0, min(top["C"], 0.24) + 1e-9, 0.01), 3))
        GRID_UQ = dict(GRID_MAIN)
    if args.production:
        OUT = os.path.join(HERE, "production")
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    log = open(os.path.join(OUT, "run.log"), "a", buffering=1)

    def say(msg):
        line = f"[{time.time() - t0:7.0f} s] {msg}"
        print(line, flush=True)
        log.write(line + "\n")

    say(f"caps {caps or 'none'}; start, {args.workers} workers, grids main={ {k: len(v) for k, v in GRID_MAIN.items()} } "
        f"uq={ {k: len(v) for k, v in GRID_UQ.items()} }, draws={n_draws}")
    draws = draw_params(n_draws)
    check_replay(draws)
    # Carbon's D is not sampled; draw_params replays the old fixed 6.3e-5, so use the production value.
    d_c = json.load(open(os.path.join(HERE, "partition_params.json")))["partition_coefficients"]["C"]["D"]
    for d in draws:
        d["d_c"] = d_c
    if args.anchor_sd is not None:  # production spread is 0.20 in both; rescale the same deviates
        for d in draws:
            d["liq_rho"] = RHO0 + (d["liq_rho"] - RHO0) * args.anchor_sd * RHO0 / 0.20
            d["liq_vphi"] = VPHI0 + (d["liq_vphi"] - VPHI0) * args.anchor_sd * VPHI0 / 0.20

    # Stage A + sensitivity: main grid at each scenario (B8), then B1 in its own pool.
    scen = [("baseline_5800K", 5800, None)] + \
           [(f"T_ICB_{t}K", t, None) for t in (5500, 5700, 6000, 6200)] + \
           [("D_Si_0.67", 5800, 0.67), ("D_Si_1.00", 5800, 1.00)]
    summ = []
    with Pool(args.workers, initializer=_init, initargs=(None, args.carbon, args.independent, caps)) as pool:
        for label, df in pool.imap_unordered(task_sweep, [(l, t, d, GRID_MAIN) for l, t, d in scen]):
            df.to_csv(os.path.join(OUT, f"sweep_{label}.csv"), index=False)
            summ.append(dict(scenario=label, **summarize(df)))
            say(f"sweep {label}: {summ[-1].get('n_admissible')} admissible, "
                f"oxygen_required={summ[-1].get('oxygen_required')}")
            if label == "baseline_5800K" and args.legacy_checks and not args.smoke:
                old = df[(df.liq_x_c <= 0.06 + 1e-9) & (df.liq_x_o <= 0.20 + 1e-9)]
                n_old = int(old.admissible.sum())
                say(f"regression: old search box gives {n_old} admissible (production 57)")
                if n_old != 57:
                    raise SystemExit("baseline does not reproduce the production sweep; stopping")
    with Pool(1, initializer=_init, initargs=("B1", args.carbon, args.independent, caps)) as pool:
        label, df = pool.apply(task_sweep, (("B1_oxide", 5800, None, GRID_MAIN),))
        df.to_csv(os.path.join(OUT, f"sweep_{label}.csv"), index=False)
        summ.append(dict(scenario=label, **summarize(df)))
        say(f"sweep {label}: {summ[-1].get('n_admissible')} admissible")
    pd.DataFrame(summ).to_csv(os.path.join(OUT, "sensitivity_summary.csv"), index=False)
    if args.production:
        write_best_fit(pd.read_csv(os.path.join(OUT, "sweep_baseline_5800K.csv")), caps, OUT)
        say("wrote best_fit.json")

    # Monte Carlo.
    rows = []
    with Pool(args.workers, initializer=_init, initargs=(None, args.carbon, args.independent, caps)) as pool:
        for i, row in enumerate(pool.imap_unordered(task_draw, [(p, GRID_UQ) for p in draws], chunksize=2), 1):
            rows.append(row)
            if i % 20 == 0 or i == n_draws:
                say(f"monte carlo {i}/{n_draws}")
                pd.DataFrame(rows).sort_values("draw").to_csv(os.path.join(OUT, "uq_draws.csv"), index=False)
    uq = pd.DataFrame(rows).sort_values("draw")
    uq.to_csv(os.path.join(OUT, "uq_draws.csv"), index=False)
    if args.legacy_checks and not args.smoke:
        ref = pd.read_csv(os.path.join(HERE, "icb_partition_uq_results.csv"))
        m = uq.merge(ref[["draw", "oc_o", "oc_s"]], on="draw")
        same = np.isclose(m.os_only_best_o.fillna(-1), m.oc_o.fillna(-1)) &             np.isclose(m.os_only_best_s.fillna(-1), m.oc_s.fillna(-1))
        say(f"regression: O+S-only best fit matches production in {int(same.sum())}/{len(m)} draws")
    a = uq[uq.admissible == True]
    pct = lambda s: [float(np.percentile(s, q)) for q in (16, 50, 84)]
    res = dict(n_draws=n_draws, frac_admissible=len(a) / len(uq),
               frac_oxygen_required=float(a.oxygen_required.mean()) if len(a) else None,
               best_o_pctl=pct(100 * a.best_o) if len(a) else None,
               best_c_pctl=pct(100 * a.best_c) if len(a) else None,
               best_si_pctl=pct(100 * a.best_si) if len(a) else None,
               best_s_pctl=pct(100 * a.best_s) if len(a) else None,
               sol_s_pctl=pct(100 * a.best_sol_s) if len(a) else None,
               jump_pctl=pct(a.best_jump) if len(a) else None)
    json.dump(res, open(os.path.join(OUT, "uq_summary.json"), "w"), indent=1)
    say(f"done: {json.dumps(res)}")


if __name__ == "__main__":
    main()
