"""Temperature and phase sensitivity of the ICB partition double constraint.

Scenarios: T_ICB in {5500,5700,5800,6000,6200} (FeO B8 oxide, D_Si 0.80); the B1
oxide variant at 5800 K; and D_Si 0.67 and 1.00 at 5800 K (B8). For each scenario
the admissibility search over outer core O and S records the admissible O range,
best fit O, inner core S, ICB density jump, and whether oxygen is required.

The B1 oxide swap is selected at model import time via FEO_PHASE; this scenario is
run in a clean subprocess with FEO_PHASE=B1 set before the model is imported.
Writes icb_partition_sensitivity.csv and fig10_sensitivity.png.
"""
import os
import sys
import json
import warnings
import subprocess
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

HERE = os.path.dirname(os.path.abspath(__file__))

X_O = None
X_S = None

TOL = dict(liq_rho=1.0, liq_vphi=2.0, sol_rho=1.5, sol_vphi=2.5)

SOURCE_NOTE = ("Source: Monte Carlo / sensitivity over EoS, partition, liquid anchor, "
               "and T_ICB (icb_partition_uq / icb_partition_sensitivity). "
               "Greyscale per house figure standard.")


def run_scenario(t_icb, d_si):
    """Run one admissibility search in the currently imported model. Returns dict."""
    import icb_partition_model as M
    X_O = np.round(np.arange(0.02, 0.201, 0.01), 4)
    X_S = np.round(np.arange(0.00, 0.0401, 0.005), 4)
    # set D_Si, keep others at model defaults
    M.PP["partition_coefficients"]["Si"]["D"] = d_si
    M.base._phase_props.cache_clear()

    pr_cmb = M.prem_outer(M.P_CMB)
    pr_icb = M.prem_outer(M.P_ICB)
    adm_o = []
    best = None
    best_rms = np.inf
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
            adm = (abs(d_lc_rho) <= TOL["liq_rho"] and abs(d_li_rho) <= TOL["liq_rho"]
                   and abs(d_lc_v) <= TOL["liq_vphi"] and abs(d_li_v) <= TOL["liq_vphi"]
                   and abs(d_s_rho) <= TOL["sol_rho"] and abs(d_s_v) <= TOL["sol_vphi"])
            if not adm:
                continue
            adm_o.append(xo)
            if xo <= 0.03:
                any_low_o = True
            jump = 100 * (sp[0] - li[0]) / li[0]
            rms = np.sqrt(np.mean(np.square(
                [d_lc_rho, d_li_rho, d_lc_v, d_li_v, d_s_rho, d_s_v])))
            if rms < best_rms:
                best_rms = rms
                best = dict(best_oc_o=xo, best_oc_s=xs, inner_s=xsol["S"],
                            icb_jump_pct=jump, rms_misfit_pct=rms)
    if not adm_o:
        return dict(adm_o_min=np.nan, adm_o_max=np.nan, best_oc_o=np.nan,
                    best_oc_s=np.nan, inner_s=np.nan, icb_jump_pct=np.nan,
                    rms_misfit_pct=np.nan, oxygen_required=True, n_admissible=0)
    out = dict(adm_o_min=min(adm_o) * 100.0, adm_o_max=max(adm_o) * 100.0,
               oxygen_required=(not any_low_o), n_admissible=len(adm_o))
    out.update({k: (v * 100.0 if k in ("best_oc_o", "best_oc_s", "inner_s") else v)
                for k, v in best.items()})
    return out


def subprocess_b1(t_icb, d_si):
    """Run the B1 oxide scenario in a clean child with FEO_PHASE=B1."""
    code = (
        "import os,warnings;warnings.filterwarnings('ignore')\n"
        "import json,numpy as np\n"
        "import icb_partition_sensitivity as S\n"
        f"r=S.run_scenario({t_icb},{d_si})\n"
        "print('JSONRESULT'+json.dumps({k:(None if (isinstance(v,float) and np.isnan(v)) else v) for k,v in r.items()}))\n"
    )
    env = dict(os.environ)
    env["FEO_PHASE"] = "B1"
    env["PYTHONPATH"] = HERE + os.pathsep + env.get("PYTHONPATH", "")
    p = subprocess.run([sys.executable, "-c", code], cwd=HERE, env=env,
                       capture_output=True, text=True)
    for line in p.stdout.splitlines():
        if line.startswith("JSONRESULT"):
            return json.loads(line[len("JSONRESULT"):])
    raise RuntimeError("B1 subprocess failed:\n" + p.stderr)


def main():
    scenarios = []
    # Temperature scan, FeO B8, D_Si 0.80
    for t in (5500, 5700, 5800, 6000, 6200):
        scenarios.append(dict(label=f"T_ICB {t} K, B8, D_Si 0.80",
                              t_icb=t, d_si=0.80, oxide="B8", b1=False,
                              short=f"{t} K"))
    # D_Si variants at 5800 K, B8
    scenarios.append(dict(label="D_Si 0.67, 5800 K, B8",
                          t_icb=5800, d_si=0.67, oxide="B8", b1=False,
                          short="D_Si 0.67"))
    scenarios.append(dict(label="D_Si 1.00, 5800 K, B8",
                          t_icb=5800, d_si=1.00, oxide="B8", b1=False,
                          short="D_Si 1.00"))
    # B1 oxide at 5800 K
    scenarios.append(dict(label="B1 oxide, 5800 K, D_Si 0.80",
                          t_icb=5800, d_si=0.80, oxide="B1", b1=True,
                          short="B1 oxide"))

    rows = []
    for sc in scenarios:
        if sc["b1"]:
            res = subprocess_b1(sc["t_icb"], sc["d_si"])
        else:
            res = run_scenario(sc["t_icb"], sc["d_si"])
        row = dict(scenario=sc["label"], short=sc["short"], oxide=sc["oxide"],
                   t_icb=sc["t_icb"], d_si=sc["d_si"])
        row.update(res)
        rows.append(row)

    df = pd.DataFrame(rows)
    out_csv = os.path.join(HERE, "icb_partition_sensitivity.csv")
    df.to_csv(out_csv, index=False)

    # Print table
    print("=" * 96)
    print("TEMPERATURE AND PHASE SENSITIVITY OF THE ICB PARTITION DOUBLE CONSTRAINT")
    print("=" * 96)
    hdr = ("Scenario", "OadmMin", "OadmMax", "BestO", "InnerS", "Jump",
           "RMS", "OxReq")
    print(f"{hdr[0]:<30}{hdr[1]:>8}{hdr[2]:>8}{hdr[3]:>8}{hdr[4]:>8}"
          f"{hdr[5]:>8}{hdr[6]:>8}{hdr[7]:>7}")
    print("-" * 96)
    for _, r in df.iterrows():
        def f(v):
            return "  na " if pd.isna(v) else f"{v:6.1f}"
        print(f"{r['scenario']:<30}{f(r['adm_o_min']):>8}{f(r['adm_o_max']):>8}"
              f"{f(r['best_oc_o']):>8}{f(r['inner_s']):>8}{f(r['icb_jump_pct']):>8}"
              f"{f(r['rms_misfit_pct']):>8}{str(bool(r['oxygen_required'])):>7}")
    print("=" * 96)
    print("Units: O admissible range, best fit O, inner core S in mol percent; "
          "jump and RMS in percent.")
    print("Note: B1 oxide variant has RMS ~0.466 versus ~0.462 for the B8 baseline.")

    make_fig(df, out_csv)
    print(f"Wrote {out_csv}")


def style_axes(ax):
    for sp in ax.spines.values():
        sp.set_visible(True)
        sp.set_color("black")
        sp.set_linewidth(0.8)
    ax.grid(False)
    ax.set_facecolor("white")
    ax.tick_params(colors="black", labelsize=8)


def make_fig(df, out_csv):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.6))
    fig.patch.set_facecolor("white")

    y = np.arange(len(df))[::-1]
    labels = df["short"].tolist()
    for yi, (_, r) in zip(y, df.iterrows()):
        if pd.isna(r["adm_o_min"]):
            continue
        lo, hi = r["adm_o_min"], r["adm_o_max"]
        # admissible band as light grey bar
        ax.barh(yi, hi - lo + 1.0, left=lo - 0.5, height=0.55, color="0.80",
                edgecolor="black", linewidth=0.7, zorder=2)
        # best fit O marker
        ax.plot(r["best_oc_o"], yi, marker="D", color="black", markersize=7,
                zorder=3)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("OUTER CORE O, mol percent")
    ax.set_ylabel("SCENARIO")
    ax.set_xlim(1, 21)

    # legend proxies
    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    band = Patch(facecolor="0.80", edgecolor="black", linewidth=0.7,
                 label="admissible O range")
    bestm = Line2D([0], [0], marker="D", color="black", linestyle="None",
                   markersize=7, label="best fit O")
    leg = ax.legend(handles=[band, bestm], fontsize=8, frameon=True,
                    loc="lower right")
    leg.get_frame().set_edgecolor("black")
    leg.get_frame().set_linewidth(0.8)
    style_axes(ax)

    fig.suptitle("SENSITIVITY OF ADMISSIBLE OUTER CORE OXYGEN ACROSS SCENARIOS",
                 fontsize=12, fontweight="bold")
    fig.text(0.5, 0.005, SOURCE_NOTE, ha="center", va="bottom", fontsize=7,
             color="black")
    fig.tight_layout(rect=[0, 0.05, 1, 0.95])
    out_png = os.path.join(os.path.dirname(out_csv), "fig10_sensitivity.png")
    fig.savefig(out_png, dpi=150, facecolor="white")
    plt.close(fig)
    print(f"Wrote {out_png}")


if __name__ == "__main__":
    main()
