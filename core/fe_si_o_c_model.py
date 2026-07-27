"""
Fe-Si-O-C inner-core QUATERNARY solid property model.

Builds BurnMan Mie-Grueneisen-Debye (mgd3) minerals from the sourced ledger
(eos_params.json), allocates a bulk Fe-Si-O-C composition into structural phases
by mass balance, and computes aggregate density, K_S, and bulk sound velocity
(Vphi) at inner-core P-T by Voigt-Reuss-Hill averaging. Density and Vphi need no
shear modulus and are computed honestly for the aggregate.

Run:  python fe_si_o_c_model.py
"""
import json
import os
from functools import lru_cache

import numpy as np
from burnman import Mineral

HERE = os.path.dirname(os.path.abspath(__file__))
FEO_PHASE = os.environ.get("FEO_PHASE", "B8").upper()
USE_FE_THERMAL_CORRECTION = os.environ.get("FE_THERMAL", "1") != "0"


def load_ledger():
    with open(os.path.join(HERE, "eos_params.json")) as f:
        return json.load(f)


def _mineral(em):
    p = dict(
        name=em["structure"], equation_of_state="mgd3",
        V_0=em["V_0"]["value"] * 1e-6, K_0=em["K_0"]["value"] * 1e9,
        Kprime_0=em["Kprime_0"]["value"], Debye_0=em["Debye_0"]["value"],
        grueneisen_0=em["grueneisen_0"]["value"], q_0=em["q_0"]["value"],
        molar_mass=em["molar_mass"]["value"], n=em["molar_mass"]["n_atoms"],
        T_0=300.0, F_0=0.0)
    sh = em.get("shear", {})
    if sh.get("G_0") is not None:
        p["G_0"] = sh["G_0"] * 1e9
        p["Gprime_0"] = sh.get("Gprime_0", 0.0)
    return Mineral(params=p)


def build_phases():
    em = load_ledger()["endmembers"]
    feo_key = "B1_FeO" if FEO_PHASE == "B1" else "B8_FeO"
    return {"hcp_fe": _mineral(em["hcp_Fe"]), "fesi": _mineral(em["B2_FeSi"]),
            "feo": _mineral(em[feo_key]), "fe7c3": _mineral(em["Fe7C3"])}


def phase_amounts(x_c, x_si, x_o):
    n_hcp = 1.0 - (10.0 / 3.0) * x_c - 2.0 * x_si - 2.0 * x_o
    if n_hcp < 0.0:
        return None
    amt = {"hcp_fe": n_hcp, "fe7c3": x_c / 3.0, "fesi": x_si, "feo": x_o}
    return amt if sum(amt.values()) > 0 else None


@lru_cache(maxsize=8192)
def _phase_props(phase_key, p_pa, t_k):
    m = _PHASES[phase_key]
    if phase_key == "hcp_fe" and USE_FE_THERMAL_CORRECTION:
        rho, ks_gpa, _ = fe_corrected_props(p_pa, t_k, m)
        return rho, ks_gpa * 1e9, m.params["molar_mass"] / rho
    m.set_state(p_pa, t_k)
    return m.density, m.isentropic_bulk_modulus_reuss, m.molar_volume


def aggregate_props(amounts, p_pa, t_k):
    mass = vol = 0.0
    parts = []
    for key, n in amounts.items():
        if n <= 0.0:
            continue
        rho_i, k_i, v_i = _phase_props(key, p_pa, t_k)
        mass += n * _PHASES[key].params["molar_mass"]
        vol += n * v_i
        parts.append((n * v_i, k_i))
    rho = mass / vol
    phi = np.array([v for v, _ in parts]) / vol
    k = np.array([kk for _, kk in parts])
    k_vrh = 0.5 * (float(np.sum(phi * k)) + float(1.0 / np.sum(phi / k)))
    return rho, k_vrh, np.sqrt(k_vrh / rho) / 1000.0


def prem_inner_core(P_gpa):
    x2 = (363.85 - np.asarray(P_gpa, dtype=float)) / 952.0
    x2 = np.clip(x2, 0.0, 0.036760)
    rho = 13.0885 - 8.8381 * x2
    vp = 11.2622 - 6.3640 * x2
    vs = 3.6678 - 4.4475 * x2
    vphi = np.sqrt(vp**2 - (4.0 / 3.0) * vs**2)
    return rho, vphi, vp, vs


_DEW = dict(V0=6.753e-6, a0=3.7e-5, m=1.87, e0=1.95e-4, g=1.339)


def _fe_extra_pressure(V, T):
    R = 8.314462
    x = V / _DEW["V0"]
    P_anh = 1.5 * R * _DEW["a0"] * _DEW["m"] * T**2 * x**_DEW["m"] / V
    P_el = 1.5 * R * _DEW["e0"] * _DEW["g"] * T**2 * x**_DEW["g"] / V
    return P_anh + P_el


def fe_corrected_props(P_pa, T_k, fe=None):
    from scipy.optimize import brentq
    if fe is None:
        fe = _PHASES["hcp_fe"]

    def Ptot(V):
        return fe.method.pressure(T_k, V, fe.params) + _fe_extra_pressure(V, T_k)

    V = brentq(lambda v: Ptot(v) - P_pa, 3.4e-6, 6.9e-6)
    rho = fe.params["molar_mass"] / V
    dV = V * 1e-4
    K_T = -V * (Ptot(V + dV) - Ptot(V - dV)) / (2 * dV)
    try:
        fe.set_state_with_volume(V, T_k)
    except Exception:
        fe.set_state(P_pa, T_k)
    K_S = K_T * (1.0 + fe.alpha * fe.grueneisen_parameter * T_k)
    return rho, K_S / 1e9, np.sqrt(K_S / rho) / 1000.0


_PHASES = build_phases()

if __name__ == "__main__":
    P_icb = 328.85e9
    for tag, (xc, xsi, xo) in {"pure Fe": (0, 0, 0), "5% O": (0, 0, 0.05)}.items():
        amt = phase_amounts(xc, xsi, xo)
        rho, ks, vphi = aggregate_props(amt, P_icb, 6000.0)
        prho, pvphi, _, _ = prem_inner_core(328.85)
        print(f"{tag}: rho={rho/1000:.3f} ({100*(rho/1000-prho)/prho:+.2f}%) Vphi={vphi:.2f}")
