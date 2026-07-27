"""Continuous (multi-stage) core formation on the sourced Fischer/Frost partitioning.

The single-stage ledger (oxygen_solver) asks what one effective equilibration is
required. This extends it to a continuous history: the Earth accretes in N mass
increments, each delivering metal that equilibrates with the current magma ocean at
the effective nickel-cobalt pressure, on the FeO-depressed liquidus. Oxygen entering
the metal is drawn from mantle FeO, so the mantle oxidation state and the core oxygen
content evolve together under a single FeO mass balance. Silicon entering the metal
returns FeO to the silicate (SiO2 + 2Fe -> Si + 2FeO), a small back-reaction on the
oxidized path, carried for rigor.

The free input is the redox history of the accreting material, c_in(phi). Two end
members bracket the literature: oxidized-early (Siebert 2013, mantle FeO high early,
drawn down as oxygen is swept to the core) and reduced-early (Wade and Wood 2005;
Rubie et al. 2011, accreting material reduced first, oxidized late).

Questions the single stage cannot answer:
  - does a 5 wt% O core emerge with a mantle FeO concentration that never exceeds the
    independently constrained early value (Siebert 2013, X_FeO ~ 0.21 ~ 28 wt%), so the
    single-stage "28 wt% local vs 19 wt% whole-mantle" tension is only a mass-average?
  - which oxidation history does a 5 wt% O core select?
"""
import numpy as np

import oxygen_solver as S

PB_EFF = 54.0                                            # effective Ni/Co pressure, Fischer 2015
FEO_NOW = S.P["mantle_FeO"]["present_wt_pct"]            # 8.05 wt%
X_CORE = S.P["mantle_magma_prediction"]["core_mass_frac"]       # 0.325
X_MANT = S.P["mantle_magma_prediction"]["mantle_mass_frac"]     # 0.675
M = S.P["molar_mass_g_mol"]
FE_PER_O = M["FeO"] / M["O"]                             # FeO mass per unit O mass, ~4.49
FEO_PER_SI = 2.0 * M["FeO"] / M["Si"]                    # FeO mass returned per unit Si mass


def _mantle_feo_wt(feo_mass, mant_mass):
    """FeO concentration (wt%) from FeO mass and total mantle mass (arbitrary units)."""
    return 100.0 * feo_mass / mant_mass if mant_mass > 0 else 0.0


def run(c_in, n=600, feo_seed=None, equilibrate=True):
    """Integrate one accretion history.

    c_in: callable phi in [0,1] -> incoming silicate FeO concentration (wt%), the redox
          history of accreting material.
    Returns per-step arrays and the final core composition.

    Bookkeeping is in mass units of one Earth mass (M_E = 1). Each of n equal increments
    adds dM = 1/n of planet mass, split X_CORE metal / X_MANT silicate.
    """
    dM = 1.0 / n
    dMc = X_CORE * dM                                   # metal to core per step
    dMs = X_MANT * dM                                   # silicate to mantle per step
    feo_seed = c_in(0.0) if feo_seed is None else feo_seed

    mant_mass = 0.0
    feo_mass = 0.0
    core_mass = 0.0
    core_O = 0.0
    core_Si = 0.0

    phis = np.linspace(1.0 / n, 1.0, n)
    a_feo, a_wO, a_wSi, a_T = (np.zeros(n) for _ in range(4))

    for i, phi in enumerate(phis):
        # incoming silicate joins the mantle first, carrying FeO at the accreting redox
        mant_mass += dMs
        feo_mass += c_in(phi) / 100.0 * dMs              # FeO mass = (wt%/100) * silicate mass
        feo_wt = _mantle_feo_wt(feo_mass, mant_mass)
        feo_wt = min(feo_wt, 45.0)                       # clamp to solver-safe range
        T = S.liquidus_T(feo_wt)
        if equilibrate:
            m = S.solve_metal(T, PB_EFF, feo_wt, mode="table2")
            wO, wSi = m["wt"]["O"], m["wt"]["Si"]
        else:
            wO = wSi = 0.0

        # oxygen to the core is drawn from mantle FeO; silicon returns FeO to the mantle
        o_to_core = wO / 100.0 * dMc
        si_to_core = wSi / 100.0 * dMc
        feo_removed = o_to_core * FE_PER_O
        feo_added_by_si = si_to_core * FEO_PER_SI
        feo_mass += (-feo_removed + feo_added_by_si)
        feo_mass = max(feo_mass, 0.0)

        core_mass += dMc
        core_O += o_to_core
        core_Si += si_to_core

        a_feo[i], a_wO[i], a_wSi[i], a_T[i] = feo_wt, wO, wSi, T

    return dict(
        phis=phis, feo=a_feo, wO=a_wO, wSi=a_wSi, T=a_T,
        core_O=100.0 * core_O / core_mass,
        core_Si=100.0 * core_Si / core_mass,
        feo_final=a_feo[-1], feo_peak=float(a_feo.max()),
        feo_removed_bulk=100.0 * (core_O * FE_PER_O) / mant_mass,   # O-equiv FeO over present mantle
    )


def const_path(feo):
    return lambda phi: feo


def oxidized_path(feo_early, feo_now=FEO_NOW, shape=1.0):
    """Accreting material oxidized early, declining to present. shape>1 front-loads it."""
    return lambda phi: feo_now + (feo_early - feo_now) * (1.0 - phi) ** shape


def reduced_path(feo_start, feo_late, shape=1.0):
    """Reduced-early (Rubie / Wade and Wood): rises from a low start to an oxidized late value."""
    return lambda phi: feo_start + (feo_late - feo_start) * phi ** shape


def solve_oxidized_for_core_O(target=5.0, shape=1.0, n=400):
    """Early-mantle FeO (wt%) whose oxidized-early history yields the target core O."""
    from scipy.optimize import brentq

    def f(fe):
        return run(oxidized_path(fe, shape=shape), n=n)["core_O"] - target
    return brentq(f, FEO_NOW + 0.5, 55.0)


if __name__ == "__main__":
    # The robust use of this model in the manuscript is directional, not a fitted number:
    # an oxidized-early accretion history builds an oxygen-rich, silicon-poor core, while a
    # reduced-early history drives silicon into the metal instead. The exact O/Si split is
    # sensitive to the accretion and equilibration profile and is NOT quoted in the paper.
    print(f"present-day mantle FeO {FEO_NOW} wt%,  core mass fraction {X_CORE}\n")

    ox = run(oxidized_path(30.0), n=600)
    red = run(reduced_path(0.5, FEO_NOW), n=600)
    print("history         core O (wt%)   core Si (wt%)   light element")
    print(f"oxidized-early   {ox['core_O']:6.2f}        {ox['core_Si']:6.2f}       oxygen-rich, Si-poor")
    print(f"reduced-early    {red['core_O']:6.2f}        {red['core_Si']:6.2f}       silicon-dominated")
    print("\nDirection is set by the partitioning (high mantle FeO suppresses Si, promotes O);")
    print("magnitude is profile-dependent, so only the direction enters the manuscript.")
