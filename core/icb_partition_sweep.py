"""Double-constraint composition search (Stage 2). Writes icb_partition_results[_carbide].csv."""
import os, sys, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import icb_partition_model as M
HERE = os.path.dirname(os.path.abspath(__file__))
X_O = np.round(np.arange(0.00, 0.205, 0.01), 3)
X_SI = np.round(np.arange(0.00, 0.105, 0.01), 3)
X_S = np.round(np.arange(0.00, 0.0605, 0.005), 4)
X_C = np.round(np.arange(0.00, 0.065, 0.02), 3)
TOL_LIQ_RHO = 1.0
TOL_LIQ_VPHI = 2.0  # CMB Vphi anchored to sourced liquid Fe (Kuwayama 2020)
TOL_SOL_RHO = 1.5
TOL_SOL_VPHI = 2.5


def run(scenario="default"):
    pr_cmb = M.prem_outer(M.P_CMB); pr_icb = M.prem_outer(M.P_ICB); rows = []; skipped = 0
    for xo in X_O:
        for xsi in X_SI:
            for xs in X_S:
                for xc in X_C:
                    if xo+xsi+xs+xc >= 0.6:
                        continue
                    lc = M.liquid_props(xc, xo, xsi, xs, M.P_CMB); li = M.liquid_props(xc, xo, xsi, xs, M.P_ICB)
                    if lc is None or li is None:
                        skipped += 1; continue
                    xsol = M.partition_to_solid({"C": xc, "O": xo, "Si": xsi, "S": xs}, scenario)
                    sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"], M.P_ICB*1e9, M.T_ICB)
                    if sp is None:
                        skipped += 1; continue
                    d_lc_rho = 100*(lc[0]-pr_cmb[0])/pr_cmb[0]; d_lc_v = 100*(lc[1]-pr_cmb[1])/pr_cmb[1]
                    d_li_rho = 100*(li[0]-pr_icb[0])/pr_icb[0]; d_li_v = 100*(li[1]-pr_icb[1])/pr_icb[1]
                    d_s_rho = 100*(sp[0]-M.PREM_IC["rho"])/M.PREM_IC["rho"]; d_s_v = 100*(sp[2]-M.PREM_IC["vphi"])/M.PREM_IC["vphi"]
                    jump = 100*(sp[0]-li[0])/li[0]
                    adm = (abs(d_lc_rho) <= TOL_LIQ_RHO and abs(d_li_rho) <= TOL_LIQ_RHO and abs(d_lc_v) <= TOL_LIQ_VPHI
                           and abs(d_li_v) <= TOL_LIQ_VPHI and abs(d_s_rho) <= TOL_SOL_RHO and abs(d_s_v) <= TOL_SOL_VPHI)
                    score = np.sqrt(np.mean(np.square([d_lc_rho, d_li_rho, d_lc_v, d_li_v, d_s_rho, d_s_v])))
                    rows.append(dict(liq_x_c=xc, liq_x_o=xo, liq_x_si=xsi, liq_x_s=xs, sol_x_c=xsol["C"], sol_x_o=xsol["O"], sol_x_si=xsol["Si"], sol_x_s=xsol["S"],
                        liq_rho_cmb=lc[0], liq_vphi_cmb=lc[1], liq_rho_icb=li[0], liq_vphi_icb=li[1], sol_rho_icb=sp[0], sol_vphi_icb=sp[2],
                        d_liq_rho_cmb=d_lc_rho, d_liq_vphi_cmb=d_lc_v, d_liq_rho_icb=d_li_rho, d_liq_vphi_icb=d_li_v, d_sol_rho_icb=d_s_rho, d_sol_vphi_icb=d_s_v,
                        icb_density_jump_pct=jump, rms_misfit_pct=score, admissible=adm))
    df = pd.DataFrame(rows); out = os.path.join(HERE, f"icb_partition_results{'_carbide' if scenario=='carbide' else ''}.csv")
    df.to_csv(out, index=False)
    print(f"Stage 2 ({scenario}): candidates={len(df)} admissible={int(df.admissible.sum())}")
    return df

if __name__ == "__main__":
    run("carbide" if len(sys.argv) > 1 and sys.argv[1] == "carbide" else "default")
