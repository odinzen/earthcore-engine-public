"""Circularity audit of the oxygen requirement (referee response).
Test A: can ANY oxygen-free liquid satisfy the double constraint? If none, the
oxygen requirement is structural (set by PREM), not inherited from the Badro
calibration. Test B: degrade the oxygen component off its Badro fit; is oxygen
still required? Run: python icb_circularity.py"""
import copy, numpy as np, icb_partition_model as M
TOL = dict(lc_rho=1.0, li_rho=1.0, lc_v=2.0, li_v=2.0, s_rho=1.5, s_v=2.5)
pr_cmb = M.prem_outer(M.P_CMB); pr_icb = M.prem_outer(M.P_ICB)


def misfits(xc, xo, xsi, xs):
    lc = M.liquid_props(xc, xo, xsi, xs, M.P_CMB); li = M.liquid_props(xc, xo, xsi, xs, M.P_ICB)
    if lc is None or li is None:
        return None
    xsol = M.partition_to_solid({"C": xc, "O": xo, "Si": xsi, "S": xs})
    sp = M.solid_props(xsol["C"], xsol["Si"], xsol["O"], xsol["S"], M.P_ICB*1e9, M.T_ICB)
    d = dict(lc_rho=100*(lc[0]-pr_cmb[0])/pr_cmb[0], lc_v=100*(lc[1]-pr_cmb[1])/pr_cmb[1],
             li_rho=100*(li[0]-pr_icb[0])/pr_icb[0], li_v=100*(li[1]-pr_icb[1])/pr_icb[1],
             s_rho=100*(sp[0]-M.PREM_IC["rho"])/M.PREM_IC["rho"], s_v=100*(sp[2]-M.PREM_IC["vphi"])/M.PREM_IC["vphi"])
    d["adm"] = all(abs(d[k]) <= TOL[k] for k in TOL); d["rms"] = np.sqrt(np.mean([d[k]**2 for k in TOL]))
    return d


def search(o_range):
    best = None; adm = []
    for xo in o_range:
        for xsi in np.arange(0, 0.121, 0.01):
            for xs in np.arange(0, 0.041, 0.005):
                for xc in (0.0, 0.02, 0.04):
                    if xo+xsi+xs+xc >= 0.5:
                        continue
                    d = misfits(xc, xo, xsi, xs)
                    if d is None:
                        continue
                    if d["adm"]:
                        adm.append((xo, xsi, xs, xc))
                    if best is None or d["rms"] < best[1]["rms"]:
                        best = ((xo, xsi, xs, xc), d)
    return best, adm

if __name__ == "__main__":
    print("=== TEST A: can an OXYGEN-FREE liquid satisfy the double constraint? ===")
    best0, adm0 = search([0.0]); c0, d0 = best0
    fails = [k for k in TOL if abs(d0[k]) > TOL[k]]
    print(f"oxygen-free admissible: {len(adm0)} ; best O-free RMS={d0['rms']:.2f}% ; FAILS {fails}")
    print(f"  best O-free: Si={c0[1]*100:.0f}% S={c0[2]*100:.1f}% C={c0[3]*100:.0f}% | sol_rho {d0['s_rho']:+.2f} liqICB_v {d0['li_v']:+.2f}")
    print("\n=== TEST B: degrade the oxygen calibration; is oxygen still required? ===")
    orig = copy.deepcopy(M.PP["liquid_outer_core"]); C1 = M.PP["liquid_outer_core"]["components_at_130GPa"]; C2 = M.PP["liquid_outer_core"]["components_at_329GPa"]
    VFe1, VFe2, KFe1, KFe2 = C1["Fe"]["V"], C2["Fe"]["V"], C1["Fe"]["K"], C2["Fe"]["K"]
    VO1, VO2 = orig["components_at_130GPa"]["O"]["V"], orig["components_at_329GPa"]["O"]["V"]
    KO1, KO2 = orig["components_at_130GPa"]["O"]["K"], orig["components_at_329GPa"]["O"]["K"]
    for scale in (1.0, 0.7, 0.5, 0.3, 0.0):
        C1["O"]["V"] = VFe1-(VFe1-VO1)*scale; C2["O"]["V"] = VFe2-(VFe2-VO2)*scale
        C1["O"]["K"] = KFe1-(KFe1-KO1)*scale; C2["O"]["K"] = KFe2-(KFe2-KO2)*scale
        best, adm = search(np.arange(0.0, 0.25, 0.01)); o_adm = [a[0] for a in adm]
        o_req = len(adm) > 0 and len([a for a in adm if a[0] <= 0.03]) == 0
        rng = f"{min(o_adm)*100:.0f}-{max(o_adm)*100:.0f}" if o_adm else "none"
        print(f"  O potency scale={scale:.1f}: admissible={len(adm):3d} O_range={rng} mol% oxygen_required={o_req}")
    M.PP["liquid_outer_core"] = orig
