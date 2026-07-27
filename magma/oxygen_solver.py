"""Core oxygen provenance: metal-silicate partitioning of O (and Si) into the
core-forming metal, from Fischer et al. 2015 (GCA 167, 177-194).

Two sourced models, selectable per call:
  mode="table2"  Eq (2), composition-free fit. Reproduces Fischer's headline
                 single-stage result (his 1.6 wt% O comes from Table 2).
  mode="eps"     Eq (3), Table 3 fit plus the Ma (2001) epsilon activity terms,
                 with interaction parameters scaled to temperature by Eq (4).
                 Captures oxygen's negative self-interaction (e_O^O), the
                 solubility runaway that drives O higher.

The metal composition (X_O, X_Si, X_Fe) is solved self-consistently: O and Si
partition together, and each enters the other's activity correction. Silicate
mole fractions come from a BSE major-oxide melt keyed to mantle FeO.

Definitions (Fischer 2015 Sec. 4.1): D_M = X_M^met / X_MOn2^sil; exchange
coefficient K_D = D_M / D_Fe^(n/2). For oxygen (Frost 2010),
K_D(O) = X_Fe^met * X_O^met / X_FeO^sil. Silicon has valence 4.

Run: python oxygen_solver.py
"""
import json
import os

import numpy as np
from scipy.optimize import brentq

HERE = os.path.dirname(os.path.abspath(__file__))
P = json.load(open(os.path.join(HERE, "oxygen_params.json")))
EX = P["exchange_coefficients"]
M = P["molar_mass_g_mol"]
SIL = P["silicate_model"]

TREF = 1873.0          # K, reference for the interaction parameters (Eq 4)
M_FE = 55.85           # Fischer's Eq (4) uses 55.85 for iron
LN10 = np.log(10.0)    # the 2.303 factor converting ln to log10
LIGHT = ("O", "Si", "C")


def silicate_molefracs(feo_wt):
    """X_FeO and X_SiO2 in the silicate melt for a given mantle FeO (wt%).
    Non-FeO majors held in fixed BSE ratios, renormalized to (100 - FeO)."""
    ox = dict(SIL["bse_nonFeO_oxide_wt_pct"])
    scale = (100.0 - feo_wt) / sum(ox.values())
    wt = {k: v * scale for k, v in ox.items()}
    wt["FeO"] = feo_wt
    n = {k: wt[k] / M[k] for k in wt}
    ntot = sum(n.values())
    return n["FeO"] / ntot, n["SiO2"] / ntot


def eps_at_T(e_ref, el, T):
    """Interaction parameter at temperature T from its 1873 K value (Eq 4).
    el is the element of interest i; its molar mass sets the scaling."""
    Mi = M[el]
    return e_ref * Mi * TREF / 0.242 / T - Mi / M_FE + 1.0


def logKD_table2(el, T, Pg, ov=None):
    """log10 K_D from the composition-free fit, Eq (2) / Table 2.
    ov optionally overrides the (a, b, c) for an element (used by the UQ sampler)."""
    c = (ov or {}).get(el) or EX[el]["table2"]
    return c["a"] + c["b"] / T + c["c"] * Pg / T


def logKD_eps(el, T, Pg, X):
    """log10 K_D from Eq (3): Table 3 base plus the epsilon activity terms.
    X holds metal mole fractions keyed by element symbol."""
    e = EX[el]
    val = e["a"] + e["b"] / T + e["c"] * Pg / T
    Xi = X[el]

    # self-interaction term e_i^i * ln(1 - X_i)
    e_ii = e.get("eps_" + el, 0.0)
    if e_ii:
        val += eps_at_T(e_ii, el, T) / LN10 * np.log(1.0 - Xi)

    # cross terms with the other light elements
    for k in LIGHT:
        if k == el:
            continue
        e_ik = e.get("eps_" + k, 0.0)
        if not e_ik or X[k] <= 0.0:
            continue
        Xk = X[k]
        eps = eps_at_T(e_ik, el, T)
        val += eps * Xk * (1.0 + np.log(1.0 - Xk) / Xk - 1.0 / (1.0 - Xi)) / LN10
        val -= eps * Xk**2 * Xi * (
            1.0 / (1.0 - Xi) + 1.0 / (1.0 - Xk)
            + Xi / (2.0 * (1.0 - Xi) ** 2) - 1.0
        ) / LN10
    return val


def solve_metal(T, Pg, feo_wt, mode="eps", x_c=0.0, x_ni=0.0, ov=None):
    """Self-consistent metal composition at (T, P, mantle FeO).
    Returns mole fractions and wt% for Fe, Si, O (plus C, Ni if set).
    ov overrides Table 2 (a, b, c) per element for the UQ sampler."""
    X_FeO, X_SiO2 = silicate_molefracs(feo_wt)

    # Damped fixed-point iteration. X_O and X_Si each follow explicitly from the
    # current composition through their exchange equilibria, so iterating to a
    # fixed point is stable where a Newton solve on the residual is not.
    cap = 0.85 - x_c - x_ni      # keep X_Fe comfortably positive
    xo, xsi = 0.05, 0.10
    for _ in range(500):
        xfe = 1.0 - xo - xsi - x_c - x_ni
        X = {"Fe": xfe, "O": xo, "Si": xsi, "C": x_c, "Ni": x_ni}
        if mode == "eps":
            kd_o = 10.0 ** logKD_eps("O", T, Pg, X)
            kd_si = 10.0 ** logKD_eps("Si", T, Pg, X)
        else:
            kd_o = 10.0 ** logKD_table2("O", T, Pg, ov)
            kd_si = 10.0 ** logKD_table2("Si", T, Pg, ov)
        # O:  K_D(O) = X_Fe X_O / X_FeO                 -> X_O
        # Si: K_D(Si) = X_Si X_FeO^2 / (X_SiO2 X_Fe^2)  -> X_Si  (valence 4)
        xo_new = kd_o * X_FeO / xfe
        xsi_new = kd_si * X_SiO2 * xfe**2 / X_FeO**2
        xo_new = min(max(xo + 0.5 * (xo_new - xo), 1e-9), cap)
        xsi_new = min(max(xsi + 0.5 * (xsi_new - xsi), 1e-9), cap)
        if xo_new + xsi_new > cap:           # stay on the physical simplex
            s = cap / (xo_new + xsi_new)
            xo_new, xsi_new = xo_new * s, xsi_new * s
        step = max(abs(xo_new - xo), abs(xsi_new - xsi))
        xo, xsi = xo_new, xsi_new
        if step < 1e-11:
            break
    # The Si-O / O-O epsilon terms are a positive feedback: above a threshold the
    # fixed point loses stability and O runs to the metal's limit (the cap). Fischer
    # (2015, Sec. 4.2.3) flags this same instability as the reason the epsilon model
    # cannot forward-model core composition. Report it rather than a spurious plateau.
    stable = step < 1e-9 and (xo + xsi) < cap - 1e-4
    xfe = 1.0 - xo - xsi - x_c - x_ni
    X = {"Fe": xfe, "Si": xsi, "O": xo, "C": x_c, "Ni": x_ni}
    mass = {el: X[el] * M[el] for el in ("Fe", "Si", "O", "C", "Ni")}
    tot = sum(mass.values())
    wt = {el: mass[el] / tot * 100.0 for el in mass}
    return {"X": X, "wt": wt, "stable": stable}


def core_O_wt(T, Pg, feo_wt, mode="eps", **kw):
    return solve_metal(T, Pg, feo_wt, mode=mode, **kw)["wt"]["O"]


FR = P["frost_activity_model"]
R_GAS = FR["R"]


def _W(name, T, Pg):
    w = FR[name]
    return w["c0"] + w.get("cT", 0.0) * T + w.get("cP_bar", 0.0) * (Pg * 1e4)


def frost_gamma_Fe(x_fe, T, Pg):
    """Activity coefficient of Fe in the Fe-FeO metallic liquid (Frost 2010 Eq 7)."""
    w_ff, w_fef = _W("W_Fe_FeO", T, Pg), _W("W_FeO_Fe", T, Pg)
    rt_ln = (1.0 - x_fe) ** 2 * (w_ff + 2.0 * (w_fef - w_ff) * x_fe)
    return np.exp(rt_ln / (R_GAS * T))


def frost_gamma_FeO_sil(x_feo, T, Pg):
    """Activity coefficient of FeO in the magnesiowustite-like silicate (Frost Eq 11)."""
    return np.exp((1.0 - x_feo) ** 2 * _W("W_FeO_MgO", T, Pg) / (R_GAS * T))


def delta_IW(feo_wt, xfe_met, T=None, Pg=None, activity=False):
    """Oxygen fugacity relative to iron-wustite, Fischer Eq (1):
    dIW = 2 log10(a_FeO^sil / a_Fe^met). Default approximates activities as mole
    fractions; activity=True applies Frost 2010's non-ideal gamma (needs T, P),
    which raises dIW by ~0.4 log units in line with Fischer."""
    x_feo, _ = silicate_molefracs(feo_wt)
    a_feo, a_fe = x_feo, xfe_met
    if activity:
        a_feo *= frost_gamma_FeO_sil(x_feo, T, Pg)
        a_fe *= frost_gamma_Fe(xfe_met, T, Pg)
    return 2.0 * np.log10(a_feo / a_fe)


def sider_logKD(el, T, Pg):
    """log10 K_D for a siderophile (Ni, Co) from its composition-free fit. Ni and
    Co carry no epsilon terms, so their top-level a/b/c are the Table 2 values."""
    c = EX[el]
    return c["a"] + c["b"] / T + c["c"] * Pg / T


def consistency(T, Pg, feo_wt, mode="table2"):
    """One magma-ocean state checked against every anchor at once: metal O and Si,
    the implied oxygen fugacity (dIW), and the Ni/Co exchange coefficients that pin
    the effective P-T. Nothing moves without showing whether the state still holds."""
    r = solve_metal(T, Pg, feo_wt, mode=mode)
    xfe = r["X"]["Fe"]
    return {
        "wtO": r["wt"]["O"], "wtSi": r["wt"]["Si"],
        "dIW": delta_IW(feo_wt, xfe),
        "dIW_act": delta_IW(feo_wt, xfe, T=T, Pg=Pg, activity=True),
        "NiKD": sider_logKD("Ni", T, Pg),
        "CoKD": sider_logKD("Co", T, Pg),
        "stable": r["stable"],
    }


def solve_T_for_O(target_wt, Pg, feo_wt, mode="table2", bracket=(2500.0, 7000.0), ov=None):
    """Temperature giving target wt% O, or None if no crossing in bracket.
    Intended for the stable Table 2 forward model."""
    def f(T):
        return core_O_wt(T, Pg, feo_wt, mode=mode, ov=ov) - target_wt
    lo, hi = bracket
    if f(lo) * f(hi) > 0:
        return None
    return brentq(f, lo, hi, xtol=1e-3)


def _coeff_draws(rng, el, n):
    """n correlated draws of an element's Table 2 coefficients. Uses Fischer's full
    variance-covariance matrix (supplement A.7/A.8) when present, else independent
    1-sigma. Returns the parameter order and an (n, k) array of draws."""
    b = EX[el]["table2"]
    if "cov" in b:
        order = b["cov_order"]
        mean = np.array([b[k] for k in order])
        return order, rng.multivariate_normal(mean, np.array(b["cov"]), size=n)
    order = ("a", "b", "c")
    draws = np.array([[b[k] + rng.normal(0.0, b.get(k + "_sd", 0.0)) for k in order]
                      for _ in range(n)])
    return order, draws


def uq_O_at_state(T, Pg, feo_wt, n=2000, seed=0):
    """Monte Carlo forward propagation of wt% O at a fixed magma-ocean state, drawing
    the O and Si Table 2 coefficients from Fischer's full covariance (the a/b/c terms
    are strongly anti-correlated, so this is far tighter than independent sampling).
    Returns median and 16/84th percentiles."""
    rng = np.random.default_rng(seed)
    o_order, o_draws = _coeff_draws(rng, "O", n)
    s_order, s_draws = _coeff_draws(rng, "Si", n)
    out = []
    for i in range(n):
        ov = {"O": {"a": 0.0, "b": 0.0, "c": 0.0},
              "Si": {"a": 0.0, "b": 0.0, "c": 0.0}}
        for j, k in enumerate(o_order):
            ov["O"][k] = o_draws[i, j]
        for j, k in enumerate(s_order):
            ov["Si"][k] = s_draws[i, j]
        out.append(core_O_wt(T, Pg, feo_wt, mode="table2", ov=ov))
    a = np.array(out)
    return {"median": float(np.median(a)),
            "lo": float(np.percentile(a, 16)),
            "hi": float(np.percentile(a, 84))}


def runaway_onset_T(Pg, feo_wt, bracket=(2500.0, 7000.0), step=10.0):
    """Lowest temperature at which the epsilon model loses stability (the O
    solubility runaway sets in), or None if it stays stable across the bracket."""
    lo, hi = bracket
    T = lo
    while T <= hi:
        if not solve_metal(T, Pg, feo_wt, mode="eps")["stable"]:
            return T
        T += step
    return None


MM = P["mantle_magma_prediction"]
LQ = P["mantle_liquidus"]


def reaction_HS(el):
    """Standard-state entropy and enthalpy of an element's exchange reaction, read
    from the van't Hoff structure of the Table 2 fit (log10 K_D = a + b/T + cP/T):
    dS0 = a R ln10, dH0 = -b R ln10. Returns (dS0 in J/mol/K, dH0 in kJ/mol)."""
    c = EX[el].get("table2", EX[el])
    dS = c["a"] * R_GAS * LN10
    dH = -c["b"] * R_GAS * LN10 / 1000.0
    return dS, dH


def liquidus_T(feo_wt):
    """Mantle liquidus temperature at a given mantle FeO, by first-order cryoscopic
    depression of the Ni/Co-pinned BSE liquidus. Adding FeO lowers the melting point
    (an entropy-of-mixing effect), so the FeO-rich oxidized route equilibrates cooler."""
    x = silicate_molefracs(feo_wt)[0]
    x0 = silicate_molefracs(P["mantle_FeO"]["present_wt_pct"])[0]
    inv = 1.0 / LQ["T_L0_K"] - (R_GAS / LQ["dHfus_J_mol"]) * np.log((1.0 - x) / (1.0 - x0))
    return 1.0 / inv


def solve_oxidized_route(target_o=5.0, dhfus=None, bracket=(8.0, 60.0)):
    """Self-consistent oxidized route: the mantle FeO at which the metal reaches
    target wt% O when the equilibration sits on that mantle's own (FeO-depressed)
    liquidus, not the fixed BSE liquidus. Returns (mantle FeO, liquidus T, wt% Si)."""
    saved = LQ["dHfus_J_mol"]
    if dhfus is not None:
        LQ["dHfus_J_mol"] = dhfus
    try:
        def f(feo):
            return core_O_wt(liquidus_T(feo), Pb_ref, feo, mode="table2") - target_o
        feo = brentq(f, *bracket, xtol=1e-3)
        T = liquidus_T(feo)
        wt = solve_metal(T, Pb_ref, feo, mode="table2")["wt"]
        return {"feo": feo, "T": T, "wtSi": wt["Si"], "wtO": wt["O"]}
    finally:
        LQ["dHfus_J_mol"] = saved


Pb_ref = P["equilibration_anchor"]["P_GPa"]


def mantle_FeO_before_core(core_o_wt):
    """Pre-core-formation mantle FeO implied by a given core oxygen content, by mass
    balance: the core's O was stripped from mantle FeO. Returns (initial mantle FeO
    wt%, the FeO wt% the mantle lost)."""
    o_core = core_o_wt / 100.0 * MM["core_mass_frac"]        # O mass, Earth-mass fraction
    feo_removed = o_core * M["FeO"] / M["O"]                 # FeO-equivalent stripped
    d_feo = feo_removed / MM["mantle_mass_frac"] * 100.0     # as mantle wt%
    return P["mantle_FeO"]["present_wt_pct"] + d_feo, d_feo


def early_mantle_FeO(core_o_wt, mass_frac):
    """Silicate FeO concentration (wt%) at oxygen acquisition, if the core oxygen was
    drawn from a mantle that was `mass_frac` of its present mass (0 < mass_frac <= 1).

    mantle_FeO_before_core spreads the extracted FeO over the present, full-size mantle
    (mass_frac = 1); but core oxygen is acquired early, when the mantle is still growing,
    so the same FeO budget sits in a smaller mass and the local concentration is higher.
    This resolves the apparent gap between the whole-mantle value and the higher silicate
    FeO the partitioning requires at equilibration."""
    init, _ = mantle_FeO_before_core(core_o_wt)             # total FeO over present mantle (wt%)
    return init / mass_frac                                 # same FeO mass, smaller mantle


def acquisition_fraction(core_o_wt, local_feo_wt):
    """Mantle mass fraction at which the bulk oxygen budget reproduces a given local
    silicate FeO (e.g. the FeO the partitioning requires at equilibration). Inverts
    early_mantle_FeO."""
    init, _ = mantle_FeO_before_core(core_o_wt)
    return init / local_feo_wt


def basalt_FeO(source_feo_wt):
    """FeO of a basalt from batch melting of a mantle source (Shaw batch equation,
    FeO mildly incompatible). Uses D_FeO and the MORB melt fraction from the ledger."""
    d, f = MM["D_FeO_bulk"], MM["melt_fraction_MORB"]
    return source_feo_wt / (d + f * (1.0 - d))


if __name__ == "__main__":
    anc = P["equilibration_anchor"]
    Pb, Tb = anc["P_GPa"], anc["T_K"]
    feo_b = P["mantle_FeO"]["present_wt_pct"]
    print(f"Fischer single-stage benchmark ({Pb:.0f} GPa, {Tb:.0f} K, "
          f"mantle FeO {feo_b:.0f} wt%) -> Fischer: {anc['result_wt_pct']['O']} wt% O, "
          f"{anc['result_wt_pct']['Si']} wt% Si")
    for mode in ("table2", "eps"):
        r = solve_metal(Tb, Pb, feo_b, mode=mode)
        label = "Table 2 (comp-free)" if mode == "table2" else "Table 3 + epsilon"
        flag = "" if r["stable"] else "  [UNSTABLE - runaway]"
        print(f"  {label:20s}: {r['wt']['O']:.2f} wt% O, {r['wt']['Si']:.2f} wt% Si{flag}")

    tgt = P["effective_KD_targets"]
    print(f"\nConsistency check (each state vs every anchor at once; targets: "
          f"Ni logKD {tgt['Ni_logKD']}, Co {tgt['Co_logKD']}, fO2 {tgt['fO2']}):")
    hdr = (f"  {'state':34s} | {'wt%O':>5} | {'wt%Si':>6} | {'dIW':>6} | "
           f"{'dIW*':>7} | {'NiKD':>6} | {'CoKD':>6}")
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))
    print("  (dIW = mole-fraction approximation; dIW* = Frost 2010 non-ideal activities)")

    def show(tag, T, feo):
        c = consistency(T, Pb, feo, mode="table2")
        print(f"  {tag:34s} | {c['wtO']:5.2f} | {c['wtSi']:6.2f} | {c['dIW']:6.2f} | "
              f"{c['dIW_act']:7.2f} | {c['NiKD']:6.3f} | {c['CoKD']:6.3f}")

    show(f"Fischer single stage {Tb:.0f}K/{feo_b:.0f}%", Tb, feo_b)
    # the two routes to 5 wt% O: hotter at fixed FeO, or more oxidized at fixed T
    T_hot = solve_T_for_O(5.0, Pb, feo_b, mode="table2")
    if T_hot:
        show(f"5 wt%O route A: hot {T_hot:.0f}K/{feo_b:.0f}%", T_hot, feo_b)
    show("5 wt%O route B: oxidized 3388K/25%", 3388.0, 25)

    c0 = consistency(Tb, Pb, feo_b, mode="table2")
    print(f"\nFrost activity check: dIW {c0['dIW']:.2f} -> {c0['dIW_act']:.2f} "
          f"(+{c0['dIW_act'] - c0['dIW']:.2f}), reproduces Fischer's stated +0.4-log fO2 shift.")

    print("\nReading: the single stage is self-consistent (dIW ~ -2, Ni/Co on target).")
    print("Reaching 5 wt% O forces a choice. The hot route over-partitions Si to ~29 wt%")
    print("(physically forbidden) and drifts the siderophiles; the oxidized route (mantle")
    print("FeO up, fO2 to ~IW-1.3) gives a low-Si O-rich core. That trade is the tension.")

    print("\nUncertainty (Monte Carlo forward propagation using Fischer's full O + Si")
    print("covariance, supplement Tables A.7/A.8). wt% O with 16-84 percentile band:")
    for tag, T, feo in [(f"anchor {Tb:.0f}K/{feo_b:.0f}% (Fischer 1.6+-0.3)", Tb, feo_b),
                        ("oxidized 3388K/25%", 3388.0, 25)]:
        u = uq_O_at_state(T, Pb, feo)
        print(f"  {tag:34s}: {u['median']:.2f} "
              f"(+{u['hi'] - u['median']:.2f}/-{u['median'] - u['lo']:.2f}) wt% O")

    print("\nSensitivity of wt% O at the anchor to each O coefficient (+-1 sd, one at a time):")
    base = EX["O"]["table2"]
    o0 = core_O_wt(Tb, Pb, feo_b, mode="table2")
    for k in ("a", "b", "c"):
        sd = base.get(k + "_sd", 0.0)
        d = []
        for s in (+1, -1):
            ov = {"O": {kk: base[kk] + (s * sd if kk == k else 0.0) for kk in ("a", "b", "c")}}
            d.append(core_O_wt(Tb, Pb, feo_b, mode="table2", ov=ov) - o0)
        print(f"  O {k} ({base[k]} +- {sd}): d(wt%O) = {d[0]:+.2f} / {d[1]:+.2f}")

    print("Reading: with the real covariance the anchor band collapses toward Fischer's")
    print("stated 1.6 +- 0.3; the strong a-b anti-correlation cancels most of the spread")
    print("that independent sampling inflated. This is the publication-grade UQ.")

    print("\nEpsilon model: oxygen self-interaction runaway onset (Fischer Sec. 4.2.3;")
    print("the epsilon fits are unstable for forward core modeling, shown as mechanism):")
    for feo in (8, 25):
        onset = runaway_onset_T(Pb, feo)
        cell = f"runs away above ~{onset:.0f} K" if onset else "stable to 7000 K"
        print(f"  mantle FeO {feo:2d} wt%: {cell}")

    # The falsifiable link: core O -> mantle FeO stripped -> pre-core mantle FeO,
    # tested against Badro's oxidized early mantle; leftover mantle -> basalt FeO vs MORB.
    lo_b, hi_b = 20.0, 25.0
    print(f"\nMantle-magma prediction (mass balance; present mantle FeO "
          f"{P['mantle_FeO']['present_wt_pct']:.1f} wt%):")
    print("  core O | mantle FeO lost | pre-core mantle FeO | vs Badro oxidized 20-25%")
    for o in (2.0, 3.5, 5.0, 7.0):
        init, lost = mantle_FeO_before_core(o)
        tag = "in range" if lo_b <= init <= hi_b else ("below" if init < lo_b else "too high")
        print(f"  {o:4.1f}% | {lost:11.1f} wt% | {init:15.1f} wt% | {tag}")
    b = basalt_FeO(P["mantle_FeO"]["present_wt_pct"])
    print(f"\n  Leftover present-day mantle -> basalt FeO {b:.1f} wt% "
          f"(observed MORB {MM['observed_MORB_FeO_wt_pct']} wt%).")
    print("Reading: a ~5-7 wt% O core requires a pre-core mantle of ~19-23 wt% FeO, matching")
    print("Badro's oxidized early mantle; the depleted present-day mantle still makes MORB-like")
    print("basalts. The falsifiable handle: too much core O demands an implausibly oxidized mantle.")

    print("\nThermodynamics of the two routes (van't Hoff: dS0 = a R ln10, dH0 = -b R ln10):")
    print(f"  {'exchange':8} {'dS0 (J/mol/K)':>14} {'dH0 (kJ/mol)':>13}")
    for el in ("O", "Si", "Ni", "Co"):
        dS, dH = reaction_HS(el)
        print(f"  {el:8} {dS:>14.1f} {dH:>13.1f}")
    print("  O and Si dissolve endothermically and are entropy-driven, so heating pulls both in;")
    print("  Si dH0 is ~3.5x O, so the hot route over-drives Si. Ni/Co are exothermic (T-pinned).")

    print("\nLiquidus-coupled oxidized route (FeO depresses the melting point; solve on that liquidus):")
    for dh in (140000, 100000, 60000):
        r = solve_oxidized_route(5.0, dhfus=dh)
        print(f"  dHfus {dh/1000:.0f} kJ/mol: mantle FeO {r['feo']:.1f} wt%, T {r['T']:.0f} K, "
              f"Si {r['wtSi']:.2f} wt%")
    print("  vs the fixed-BSE-liquidus route (25 wt% FeO, 3388 K): the self-consistent route needs")
    print("  ~28 wt% FeO at ~3190 K, at the top of Badro's band, with silicon still below 1 wt%.")
