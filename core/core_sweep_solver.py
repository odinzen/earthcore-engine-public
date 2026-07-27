"""Fe-Si-O-C quaternary property sweep vs PREM (Stage 1). Writes core_sweep_results.csv."""
import os, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
import fe_si_o_c_model as model
HERE = os.path.dirname(os.path.abspath(__file__))
PRESSURES_GPA = np.arange(330, 366, 2)
TEMPERATURES_K = np.arange(5000, 7250, 250)
X_C = np.linspace(0.0, 0.05, 5); X_SI = np.linspace(0.0, 0.10, 5); X_O = np.linspace(0.0, 0.10, 5)


def run():
    rows = []; skipped = 0
    for p_gpa in PRESSURES_GPA:
        p_pa = float(p_gpa)*1e9; prem_rho, prem_vphi, _, _ = model.prem_inner_core(p_gpa)
        for t_k in TEMPERATURES_K:
            t_k = float(t_k)
            for x_c in X_C:
                for x_si in X_SI:
                    for x_o in X_O:
                        amt = model.phase_amounts(x_c, x_si, x_o)
                        if amt is None:
                            skipped += 1; continue
                        try:
                            rho, k_s, vphi = model.aggregate_props(amt, p_pa, t_k)
                        except Exception:
                            skipped += 1; continue
                        rg = rho/1000.0; nt = sum(amt.values())
                        rows.append(dict(pressure_gpa=float(p_gpa), temperature_k=t_k, bulk_x_c=x_c, bulk_x_si=x_si, bulk_x_o=x_o,
                            phase_frac_hcp_fe=amt["hcp_fe"]/nt, phase_frac_fe7c3=amt["fe7c3"]/nt, phase_frac_fesi=amt["fesi"]/nt, phase_frac_feo=amt["feo"]/nt,
                            density_g_cm3=rg, bulk_modulus_gpa=k_s/1e9, v_phi_km_s=vphi, prem_density_g_cm3=float(prem_rho), prem_v_phi_km_s=float(prem_vphi),
                            density_misfit_pct=100*(rg-prem_rho)/prem_rho, v_phi_misfit_pct=100*(vphi-prem_vphi)/prem_vphi))
    df = pd.DataFrame(rows); df.to_csv(os.path.join(HERE, "core_sweep_results.csv"), index=False)
    print(f"Stage 1 sweep: valid rows={len(df)} skipped={skipped}")
    return df

if __name__ == "__main__":
    run()
