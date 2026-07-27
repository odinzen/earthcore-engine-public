"""
Chemical-potential partitioning across the inner-core boundary (Fe Ni Si O C S).

Freezing-front closure: mu_i^liquid = mu_i^solid gives D_i = x_solid/x_liquid.
D_i calibrated to the ab initio literature (Alfe 2002 O,S,Si; Li 2019 C; D_Ni ~1).
Nickel is carried at the cosmochemical 5 at% in both reservoirs and substitutes
into the hcp metal at the iron atomic volume, raising density through mass alone.
The liquid arm is a calibrated proxy (partition_params.json), not an ab initio EoS.

Run:  python icb_partition_model.py
"""
import json
import os
import numpy as np

import fe_si_o_c_model as base

HERE = os.path.dirname(os.path.abspath(__file__))
with open(os.path.join(HERE, "partition_params.json")) as f:
    PP = json.load(f)

M_ATOM = {"Fe": 0.055845, "C": 0.012011, "O": 0.015999, "Si": 0.0280855, "S": 0.03206, "Ni": 0.058693}
T_ICB = PP["validation_targets_PREM_outer_core"]["T_ICB_K"]["value"]
P_ICB, P_CMB = 328.85, 135.75
PREM_IC = {"rho": 12.764, "vphi": 10.259}
X_NI = PP["nickel"]["core_X_Ni_atomic"]

if "fe3s" not in base._PHASES:
    base._PHASES["fe3s"] = base._mineral(PP["solid_sulfide_endmember"]["Fe3S"])


def solid_phase_amounts(x_c, x_si, x_o, x_s, x_ni=X_NI):
    n_hcp = 1.0 - (10.0/3.0)*x_c - 2.0*x_si - 2.0*x_o - 4.0*x_s - x_ni
    if n_hcp < 0.0:
        return None
    return {"hcp_fe": n_hcp, "fe7c3": x_c/3.0, "fesi": x_si, "feo": x_o, "fe3s": x_s}


def solid_props(x_c, x_si, x_o, x_s, p_pa, t_k, x_ni=X_NI):
    amt = solid_phase_amounts(x_c, x_si, x_o, x_s, x_ni)
    if amt is None:
        return None
    n_fe = amt["hcp_fe"]; metal_atoms = n_fe + x_ni
    rho_fe, k_fe, v_fe = base._phase_props("hcp_fe", p_pa, t_k)
    mass = n_fe*M_ATOM["Fe"] + x_ni*M_ATOM["Ni"]
    vol = metal_atoms*v_fe
    parts = [(metal_atoms*v_fe, k_fe)]
    for key in ("fe7c3", "fesi", "feo", "fe3s"):
        n = amt[key]
        if n <= 0.0:
            continue
        rho_i, k_i, v_i = base._phase_props(key, p_pa, t_k)
        mass += n*base._PHASES[key].params["molar_mass"]; vol += n*v_i
        parts.append((n*v_i, k_i))
    rho = mass/vol
    phi = np.array([v for v, _ in parts])/vol
    k = np.array([kk for _, kk in parts])
    k_vrh = 0.5*(float(np.sum(phi*k)) + float(1.0/np.sum(phi/k)))
    return rho/1000.0, k_vrh, np.sqrt(k_vrh/rho)/1000.0


_LIQ = PP["liquid_outer_core"]
_P1, _P2 = _LIQ["anchor_pressures_GPa"]
_C1, _C2 = _LIQ["components_at_130GPa"], _LIQ["components_at_329GPa"]


def _comp(el, p_gpa):
    f = (p_gpa - _P1)/(_P2 - _P1)
    return (_C1[el]["V"] + f*(_C2[el]["V"]-_C1[el]["V"]),
            _C1[el]["K"] + f*(_C2[el]["K"]-_C1[el]["K"]))


def liquid_props(x_c, x_o, x_si, x_s, p_gpa, x_ni=X_NI):
    x = {"C": x_c, "O": x_o, "Si": x_si, "S": x_s, "Ni": x_ni}
    x["Fe"] = 1.0 - sum(x.values())
    if x["Fe"] <= 0.0:
        return None
    v_mix = m_mix = inv_k = 0.0
    for el, xi in x.items():
        if xi <= 0.0:
            continue
        v_i, k_i = _comp(el, p_gpa)
        v_mix += xi*v_i; m_mix += xi*M_ATOM[el]; inv_k += xi/k_i
    rho = m_mix/(v_mix*1e-6)
    return rho/1000.0, np.sqrt((1.0/inv_k)*1e9/rho)/1000.0


def partition_to_solid(x_liq, scenario="default"):
    pc = PP["partition_coefficients"]; out = {}
    for el in ("C", "Si", "O", "S", "Ni"):
        D = pc[el]["D"]
        if el == "C" and scenario == "carbide":
            D = pc["C"]["D_scenario_carbide"]
        out[el] = D*x_liq.get(el, 0.0)
    return out


def prem_outer(p_gpa):
    t = PP["validation_targets_PREM_outer_core"]
    p1, p2 = t["CMB"]["P_GPa"], t["ICB_liquid_side"]["P_GPa"]
    f = (p_gpa-p1)/(p2-p1)
    return (t["CMB"]["density_g_cm3"] + f*(t["ICB_liquid_side"]["density_g_cm3"]-t["CMB"]["density_g_cm3"]),
            t["CMB"]["vphi_km_s"] + f*(t["ICB_liquid_side"]["vphi_km_s"]-t["CMB"]["vphi_km_s"]))


def wt_to_atom(w):
    n = {el: w.get(el, 0.0)/(M_ATOM[el]*1000) for el in ("C", "O", "Si", "S")}
    n["Fe"] = (1.0-sum(w.values()))/(M_ATOM["Fe"]*1000)
    tot = sum(n.values())
    return {el: n[el]/tot for el in ("C", "O", "Si", "S")}


if __name__ == "__main__":
    print(f"T_ICB={T_ICB} K; Ni={X_NI*100:.0f} at%; D_Ni={PP['partition_coefficients']['Ni']['D']}")
    for tag, (xc, xo, xsi, xs) in {"pureFe@329": (0,0,0,0), "Fe.74S.03O.23@329": (0,0.23,0,0.03)}.items():
        print(f"  {tag}: liquid", [round(v,3) for v in liquid_props(xc,xo,xsi,xs,329.0,x_ni=0.0)])
    xl = wt_to_atom({"O": 0.05, "S": 0.013})
    print("  example partition O,S ->", {k: round(v*100,3) for k,v in partition_to_solid({**xl}).items()})
