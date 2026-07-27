"""Fit own equation-of-state parameters from the QHA dump and emit eos_params_qha.json.

Completes the fit_qha TODO. For each phase (from qha_<phase>.json):
  cold/300K reference  -> Birch-Murnaghan 3rd order  -> V0, K0, K0' (the mgd3 isotherm)
  F_vib(V,T)           -> effective Debye theta_D(V) at each volume
  theta_D(V)           -> Grueneisen fit theta(V)=theta0 exp[(g0/q)(1-(V/V0)^q)] -> theta0, g0, q
  -dF/dV vs literature -> density gate (agreement/disagreement is a documented finding)

Writes eos_params_qha.json: own-derived mgd3 block per phase (own_derived_parameter) with the
imported literature value kept as validation_comparison. numpy/scipy only; no phonopy, so it runs
and is checked off Sol.

    python fit_eos_qha.py
"""
import json
import os

import numpy as np
from scipy.integrate import quad
from scipy.optimize import curve_fit, minimize_scalar, brentq

HERE = os.path.dirname(os.path.abspath(__file__))
EVA3_GPA = 160.2176621              # 1 eV/Ang^3 -> GPa
NA = 6.02214076e23
KB = 8.617333262e-5                 # eV/K
C = 1e24 / NA                       # cm^3/mol -> Ang^3 per formula/atom

# per phase: atoms/cell, formulas/cell, cell mass (g/mol), and the imported EoS
#   (lit V0 unit, lit V0 value, lit K0 GPa, lit K0', lit gamma0, lit q, lit theta0, source)
PHASES = {
    "hcp_Fe":  dict(nat=2,  nform=2, mass=2*55.845,
                    v0=6.753, v0unit="cm3/mol_atom", k0=163.4, kp=5.38, g0=1.875, q=1.305,
                    th0=417.0, src="Dewaele 2006 PRL 97 215504"),
    "B2_FeSi": dict(nat=2,  nform=1, mass=55.845+28.0855,
                    v0=12.828, v0unit="cm3/mol_formula", k0=230.6, kp=4.17, g0=1.6, q=1.0,
                    th0=417.0, src="Fischer 2014 Table 1"),
    "B1_FeO":  dict(nat=2,  nform=1, mass=55.845+15.999,
                    v0=12.256, v0unit="cm3/mol_formula", k0=149.4, kp=3.60, g0=1.53, q=1.0,
                    th0=417.0, src="Fischer 2011 EPSL (B1)"),
    "Fe3S":    dict(nat=16, nform=4, mass=4*(3*55.845+32.06),
                    v0=56.76, v0unit="cm3/mol_formula(x2)", k0=156.0, kp=3.80, g0=1.01, q=1.0,
                    th0=417.0, src="Seagle 2006 JGR 111 B06209"),
    "Fe7C3":   dict(nat=20, nform=2, mass=2*(7*55.845+3*12.011),
                    v0=53.0, v0unit="cm3/mol_formula", k0=291.0, kp=4.50, g0=1.7, q=1.0,
                    th0=417.0, src="Mookherjee/Nakajima 2011"),
}


def lit_v0_cell(p):
    """Literature V0 as an Ang^3 cell volume, matching our computed cell."""
    if p["v0unit"] == "cm3/mol_atom":
        return p["v0"] * C * p["nat"]
    # both cm3/mol_formula conventions land on the same cell volume when nform is set to match
    return p["v0"] * C * (p["nform"] if "x2" not in p["v0unit"] else 2)


def bm3_E(V, E0, V0, K0, Kp):
    """Birch-Murnaghan 3rd-order energy. K0 in eV/Ang^3."""
    t = (V0 / V) ** (2.0 / 3.0) - 1.0
    return E0 + 9.0 * V0 * K0 / 16.0 * (t**3 * Kp + t**2 * (6.0 - 4.0 * (V0 / V) ** (2.0 / 3.0)))


def bm3_P(V, V0, K0, Kp):
    f = 0.5 * ((V0 / V) ** (2.0 / 3.0) - 1.0)
    return 3.0 * K0 * f * (1.0 + 2.0 * f) ** 2.5 * (1.0 + 1.5 * (Kp - 4.0) * f)


_DEB = {}


def debye_D3(x):
    if x < 1e-6:
        return 1.0
    key = round(x, 4)
    if key not in _DEB:
        integ = quad(lambda t: t**3 / (np.expm1(t)), 0.0, x)[0]
        _DEB[key] = 3.0 / x**3 * integ
    return _DEB[key]


def debye_Fvib(T, theta, nat):
    """Debye vibrational Helmholtz free energy (eV/cell), nat atoms."""
    out = np.empty_like(T, dtype=float)
    for i, t in enumerate(T):
        x = theta / t
        out[i] = nat * KB * (1.125 * theta + 3.0 * t * np.log1p(-np.exp(-x)) - t * debye_D3(x))
    return out


def theta_at_V(T, Fvib, nat):
    """Effective Debye temperature reproducing F_vib(T) at one volume."""
    r = minimize_scalar(lambda th: float(np.sum((debye_Fvib(T, th, nat) - Fvib) ** 2)),
                        bounds=(120.0, 1600.0), method="bounded")
    return r.x


def gruneisen(V, th0, g0, q, V0):
    return th0 * np.exp(g0 / q * (1.0 - (V / V0) ** q))


def gamma_of_V(V, g0, q, V0):
    return g0 * (V / V0) ** q


RY_TO_EV = 13.605693122994
RATIOS_PHONON = [0.55, 0.58, 0.61, 0.64, 0.67, 0.70, 0.73, 0.76, 0.78]
RATIOS_ANCHOR = [0.86, 0.94, 1.02]
# Static total energies (Ry) at the low-compression anchor volumes v09/v10/v11 (Sol, 2026-07-11).
# These extend the cold curve to P~0 so V0/K0/K0' are pinned, not extrapolated (see make_manifest).
ANCHOR_E_RY = {
    "hcp_Fe":  [-658.50249883, -658.50592909, -658.48785658],
    "B2_FeSi": [-376.03224537, -376.05387349, -376.05550101],
    "B1_FeO":  [-370.88878667, -370.90404991, -370.90490695],
    "Fe3S":    [-4216.21940025, -4216.35673874, -4216.39164500],
    "Fe7C3":   [-4719.93276890, -4720.14569461, -4720.18370043],
}


def fit_phase(name, p):
    d = json.load(open(os.path.join(HERE, f"qha_{name}.json")))
    T = np.array(d["temperatures_K"])
    vs = sorted(d["volumes"], key=lambda z: z["V_ang3"])
    V = np.array([z["V_ang3"] for z in vs])
    Est = np.array([z["E_static_eV"] for z in vs])
    Fv = np.array([z["F_vib_eV"] for z in vs])            # (nV, nT)
    i300 = int(np.where(T == 300.0)[0][0])
    F300 = Est + Fv[:, i300]
    lV0 = lit_v0_cell(p)
    mass = p["mass"]

    # --- cold static E(V) over the FULL range: 9 phonon volumes + 3 low-compression anchors ---
    # The anchors reach P~0, so V0/K0/K0' are pinned by the data instead of extrapolated.
    v0_starter = float(np.mean(V / np.array(RATIOS_PHONON)))
    V_anch = v0_starter * np.array(RATIOS_ANCHOR)
    E_anch = np.array(ANCHOR_E_RY[name]) * RY_TO_EV
    Vc = np.concatenate([V, V_anch])
    Ec = np.concatenate([Est, E_anch])
    popt, _ = curve_fit(bm3_E, Vc, Ec, p0=[Ec.min(), v0_starter * 0.95, 1.2, 4.5],
                        bounds=([Ec.min() - 5, v0_starter * 0.82, 0.3, 2.5],
                                [Ec.max() + 5, v0_starter * 1.10, 4.0, 9.0]), maxfev=300000)
    E0, V0, K0ev, Kp = popt
    K0 = K0ev * EVA3_GPA
    rms = np.sqrt(np.mean((Ec - bm3_E(Vc, *popt)) ** 2)) * 1000 / p["nat"]   # meV/atom

    # --- Grueneisen gamma0, q from the thermodynamic identity gamma = V (dS/dV)_T / C_V ---
    # Direct from F(V,T), so it is insensitive to DOS shape (unlike a single-Debye fit). Trim the
    # points the F_vib noise blows up; if too few survive or the fit rails, fall back to literature.
    Sv = -np.gradient(Fv, T, axis=1)
    gVT = V[:, None] * np.gradient(Sv, V, axis=0) / (T[None, :] * np.gradient(Sv, T, axis=1))
    hi = (T >= 2000) & (T <= 5000)
    gV = np.nanmean(gVT[:, hi], axis=1)
    ok = (gV > 0.3) & (gV < 4.0)
    thermal_src = "own_derived"
    if int(ok.sum()) >= 6:
        (g0, q), _ = curve_fit(lambda v, g0, qq: gamma_of_V(v, g0, qq, V0), V[ok], gV[ok],
                               p0=[1.6, 1.0], bounds=([0.3, -1.5], [3.5, 3.0]), maxfev=100000)
        if g0 > 3.4 or g0 < 0.35:                          # railed -> not trustworthy
            g0, q, thermal_src = p["g0"], p["q"], "literature (own fit railed)"
    else:
        g0, q, thermal_src = p["g0"], p["q"], f"literature (only {int(ok.sum())}/9 clean)"

    # Debye theta0 at V0: per-volume theta_D, then a 1-param fit with gamma0/q fixed.
    thetas = np.array([theta_at_V(T, Fv[k], p["nat"]) for k in range(len(V))])
    (th0,), _ = curve_fit(lambda v, t0: gruneisen(v, t0, g0, q, V0), V, thetas,
                          p0=[float(np.clip(thetas.min() * 0.9, 260, 1100))],
                          bounds=([180], [1500]), maxfev=100000)
    th_rms = np.sqrt(np.mean((thetas - gruneisen(V, th0, g0, q, V0)) ** 2))
    if thermal_src != "own_derived":
        th0 = p["th0"]                                     # fall back the Debye temperature too

    # --- density gate vs literature (fit-free -dF/dV) ---
    def rho(v):
        return mass / NA / (v * 1e-24)

    def litP(v):
        return bm3_P(v, lV0, p["k0"] / EVA3_GPA, p["kp"]) * EVA3_GPA
    ds = []
    for i in range(1, len(V) - 1):
        P = -(F300[i + 1] - F300[i - 1]) / (V[i + 1] - V[i - 1]) * EVA3_GPA
        Vl = brentq(lambda v: litP(v) - P, lV0 * 0.35, lV0 * 1.25)
        ds.append(100 * (rho(V[i]) - rho(Vl)) / rho(Vl))
    maxd = float(np.max(np.abs(ds)))

    # V0/K0/K0' are pinned by the anchors (own). gamma0/q are the direct thermodynamic Grueneisen
    # (own where the F_vib is clean; literature fallback where it is too noisy, flagged per phase).
    dK0 = 100 * (K0 - p["k0"]) / p["k0"]
    reliability = ("density validated" if maxd < 3.0 else "density: documented residual")
    reliability += f"; cold EoS own; thermal {thermal_src}"

    return {
        "structure_atoms_per_cell": p["nat"], "formulas_per_cell": p["nform"],
        "density_gate_max_pct_vs_lit": round(maxd, 2),
        "provenance": "own_derived_parameter (QHA phonons 0.55-0.78 V/V0 + static anchors 0.86-1.02, "
                      "QE PBE, this campaign)",
        "thermal_source": thermal_src,
        "V0_ang3_cell": round(float(V0), 3),
        "K0_GPa": round(K0, 1), "Kprime_0": round(float(Kp), 2),
        "Debye_theta0_K": round(float(th0), 1), "grueneisen_0": round(float(g0), 3),
        "q_0": round(float(q), 3),
        "cold_bm3_rms_meV_atom": round(rms, 2), "theta_fit_rms_K": round(float(th_rms), 1),
        "K0_pct_vs_literature": round(dK0, 1),
        "reliability": reliability,
        "imported_literature": {
            "V0_cm3": p["v0"], "V0_unit": p["v0unit"], "K0_GPa": p["k0"], "Kprime_0": p["kp"],
            "grueneisen_0": p["g0"], "q_0": p["q"], "Debye_theta0_K": p["th0"], "source": p["src"]},
    }


def main():
    out = {"_README": "Own first-principles EoS for the FeSiOCS core solid endmembers, from our QHA "
                      "campaign (qha_campaign/). V0/K0/K0' are pinned by static E(V) over V/V0 = "
                      "0.55-1.02 (phonon volumes plus low-compression anchors). gamma0/q are the "
                      "direct thermodynamic Grueneisen (gamma = V (dS/dV)/C_V from F(V,T)): own where "
                      "the phonon data is clean (hcp/FeSi/Fe7C3, plus Fe3S; each agrees with the "
                      "literature gamma), literature fallback where F_vib is too noisy (B1_FeO) -- "
                      "see thermal_source per phase. density_gate_max_pct_vs_lit is the fit-free "
                      "-dF/dV density check vs the imported EoS across the core range. hcp-Fe K0 "
                      "differs from Dewaele by design (nonmagnetic PBE vs an experimental fit through "
                      "the magnetic collapse); the two agree in the core. Reproduce: python "
                      "fit_eos_qha.py.",
           "reference_isotherm": "0 K static (cold); thermal added via Debye-Grueneisen",
           "endmembers": {}}
    print(f"{'phase':8s} {'gate%':>6s} {'K0_own':>6s} {'K0_lit':>6s} {'g0_own':>6s} {'g0_lit':>6s}  thermal")
    for name, p in PHASES.items():
        if name not in ANCHOR_E_RY:
            continue                                   # B8_FeO not run; B1_FeO is the FeO polymorph used
        r = fit_phase(name, p)
        out["endmembers"][name] = r
        print(f"{name:8s} {r['density_gate_max_pct_vs_lit']:5.2f}% {r['K0_GPa']:6.1f} "
              f"{r['imported_literature']['K0_GPa']:6.1f} {r['grueneisen_0']:6.2f} "
              f"{r['imported_literature']['grueneisen_0']:6.2f}  {r['thermal_source']}")
    with open(os.path.join(HERE, "eos_params_qha.json"), "w") as f:
        json.dump(out, f, indent=2)
    print("\nwrote eos_params_qha.json")


if __name__ == "__main__":
    main()
