"""Aggregate Vp, Vs, Poisson of the inner-core solid vs PREM. Shear moduli imported
per phase (partition_params shear_core_GPa), flagged: K and G come from different
studies. The elastic aggregate is too stiff in shear (inner core shear anomaly);
resolution lies in premelting (Martorell 2013), not composition. Run: python icb_velocity.py"""
import os, numpy as np, icb_partition_model as M
SC = M.PP["shear_core_GPa"]; P = M.P_ICB*1e9; T = M.T_ICB


def _phase_G(key):
    if key == "fe7c3":
        m = M.base._PHASES["fe7c3"]; m.set_state(P, T); return m.shear_modulus/1e9
    return SC[key]["G"]


def vp_vs(x_c, x_si, x_o, x_s, x_ni=M.X_NI):
    amt = M.solid_phase_amounts(x_c, x_si, x_o, x_s, x_ni)
    if amt is None:
        return None
    rho_gcc, K_Pa, vphi = M.solid_props(x_c, x_si, x_o, x_s, P, T, x_ni); rho = rho_gcc*1000.0
    metal_atoms = amt["hcp_fe"]+x_ni; _, _, v_fe = M.base._phase_props("hcp_fe", P, T)
    vol = {"hcp_fe": metal_atoms*v_fe}
    for key in ("fe7c3", "fesi", "feo", "fe3s"):
        n = amt[key]
        if n > 0:
            _, _, vi = M.base._phase_props(key, P, T); vol[key] = n*vi
    Vtot = sum(vol.values()); phi = {k: vol[k]/Vtot for k in vol}; Gk = {k: _phase_G(k) for k in vol}
    Gv = sum(phi[k]*Gk[k] for k in vol); Gr = 1.0/sum(phi[k]/Gk[k] for k in vol); G = 0.5*(Gv+Gr)*1e9
    vp = np.sqrt((K_Pa+4.0/3.0*G)/rho)/1000.0; vs = np.sqrt(G/rho)/1000.0
    pois = (vp**2-2*vs**2)/(2*(vp**2-vs**2))
    return dict(rho=rho_gcc, K=K_Pa/1e9, G=G/1e9, vphi=vphi, vp=vp, vs=vs, poisson=pois)


if __name__ == "__main__":
    prem = SC["PREM_inner_core"]
    print(f"PREM IC: Vp={prem['Vp_km_s']} Vs={prem['Vs_km_s']} Poisson={prem['poisson']}")
    import csv
    rows = []
    import json
    sol = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "production", "best_fit.json")))["solid_x"]
    for tag, (c, si, o, s) in {"pure Fe Ni": (0,0,0,0), "best fit solid": (sol["C"], sol["Si"], sol["O"], sol["S"])}.items():
        r = vp_vs(c, si, o, s)
        print(f"{tag}: Vp={r['vp']:.2f}({100*(r['vp']-prem['Vp_km_s'])/prem['Vp_km_s']:+.1f}%) Vphi={r['vphi']:.2f} Vs={r['vs']:.2f}({100*(r['vs']-prem['Vs_km_s'])/prem['Vs_km_s']:+.1f}%) Poisson={r['poisson']:.3f}")
        r["case"] = tag; rows.append(r)
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "icb_velocity_results.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["case", "rho", "K", "G", "vphi", "vp", "vs", "poisson"]); w.writeheader()
        for r in rows: w.writerow({k: r[k] for k in w.fieldnames})
