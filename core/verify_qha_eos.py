"""Cross-check the FeSiOCS inner-core result on our OWN first-principles EoS.

Builds a second ledger (core/eos_params_qha.json) in which the four inner-core solid phases
(hcp-Fe, B2-FeSi, B1-FeO, Fe7C3) carry the QHA campaign's own-derived mgd3 parameters
(qha_campaign/eos_params_qha.json) instead of the imported literature values, keeping the
structural fields (molar mass, atom count, shear) unchanged. Then computes the aggregate
inner-core density at the ICB for a few oxygen contents with both ledgers and compares to PREM.

The point is not to replace the production ledger (the literature EoS is validated and the QHA
reference isotherm is 0 K static, a ~1% density offset) but to confirm the oxygen-rich conclusion
is robust to swapping in our own physics. Needs burnman (environment.yml).

    python verify_qha_eos.py
"""
import json
import os
import warnings

warnings.filterwarnings("ignore")
from burnman import Mineral  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
NA = 6.02214076e23
CELL_TO_CM3MOL = NA / 1e24
UNITS_PER_CELL = {"hcp_Fe": 2, "B2_FeSi": 1, "B1_FeO": 1, "Fe7C3": 2}   # ledger units per QHA cell
KEYS = {"hcp_fe": "hcp_Fe", "fesi": "B2_FeSi", "feo": "B1_FeO", "fe7c3": "Fe7C3"}
P_ICB, T_ICB, PREM_RHO = 328.85e9, 6000.0, 12.7636   # PREM inner core at the ICB (13.0885 is the centre)


def mineral(block):
    p = dict(name=block["structure"], equation_of_state="mgd3",
             V_0=block["V_0"]["value"] * 1e-6, K_0=block["K_0"]["value"] * 1e9,
             Kprime_0=block["Kprime_0"]["value"], Debye_0=block["Debye_0"]["value"],
             grueneisen_0=block["grueneisen_0"]["value"], q_0=block["q_0"]["value"],
             molar_mass=block["molar_mass"]["value"], n=block["molar_mass"]["n_atoms"],
             T_0=300.0, F_0=0.0)
    return Mineral(params=p)


def own_block(lit_block, name, own):
    """Literature block with the EoS parameters overwritten by our own-derived values."""
    b = json.loads(json.dumps(lit_block))
    o = own[name]
    b["V_0"] = {"value": round(o["V0_ang3_cell"] * CELL_TO_CM3MOL / UNITS_PER_CELL[name], 4),
                "unit": lit_block["V_0"]["unit"], "source": "QHA campaign (own)"}
    b["K_0"] = {"value": o["K0_GPa"], "unit": "GPa", "source": "QHA campaign (own)"}
    b["Kprime_0"] = {"value": o["Kprime_0"], "unit": "-", "source": "QHA campaign (own)"}
    b["Debye_0"] = {"value": o["Debye_theta0_K"], "unit": "K", "source": f"QHA ({o['thermal_source']})"}
    b["grueneisen_0"] = {"value": o["grueneisen_0"], "unit": "-", "source": f"QHA ({o['thermal_source']})"}
    b["q_0"] = {"value": o["q_0"], "unit": "-", "source": f"QHA ({o['thermal_source']})"}
    b["validation_comparison_literature"] = {
        "V0": lit_block["V_0"], "K0": lit_block["K_0"], "Kprime_0": lit_block["Kprime_0"],
        "grueneisen_0": lit_block["grueneisen_0"], "q_0": lit_block["q_0"]}
    return b


def aggregate_density(phases, amounts):
    mass = vol = 0.0
    for key, n in amounts.items():
        if n <= 0.0:
            continue
        m = phases[key]
        m.set_state(P_ICB, T_ICB)
        mass += n * m.params["molar_mass"]
        vol += n * m.molar_volume
    return mass / vol / 1000.0


def main():
    lit = json.load(open(os.path.join(HERE, "eos_params.json")))
    own = json.load(open(os.path.join(HERE, "..", "qha_campaign", "eos_params_qha.json")))["endmembers"]
    em = lit["endmembers"]

    qha_ledger = json.loads(json.dumps(lit))
    for name in UNITS_PER_CELL:
        qha_ledger["endmembers"][name] = own_block(em[name], name, own)
    qha_ledger["_README"] = ("Inner-core ledger with the four solid phases on the QHA own-derived "
                             "mgd3 EoS (literature kept per block as validation_comparison_literature). "
                             "Cross-check only; see verify_qha_eos.py. Production stays eos_params.json.")
    with open(os.path.join(HERE, "eos_params_qha.json"), "w") as f:
        json.dump(qha_ledger, f, indent=2)

    lit_phases = {k: mineral(em[v]) for k, v in KEYS.items()}
    own_phases = {k: mineral(qha_ledger["endmembers"][v]) for k, v in KEYS.items()}

    print(f"Inner core at the ICB ({P_ICB/1e9:.1f} GPa, {T_ICB:.0f} K), FeO = B1, MGD base.")
    print(f"PREM density at the top of the inner core: {PREM_RHO:.3f} g/cm3\n")
    print("No iron thermal correction here, so compare the two ledgers with each other, not with Stage 1.")
    print(f"{'oxygen':8s} {'rho_lit':>8s} {'rho_own':>8s} {'lit vs PREM':>12s} {'own vs PREM':>12s} {'own vs lit':>11s}")
    for xo in (0.0, 0.03, 0.05, 0.08):
        n_hcp = 1.0 - 2.0 * xo
        amt = {"hcp_fe": n_hcp, "feo": xo, "fesi": 0.0, "fe7c3": 0.0}
        rl = aggregate_density(lit_phases, amt)
        ro = aggregate_density(own_phases, amt)
        print(f"{xo*100:5.0f} %  {rl:8.3f} {ro:8.3f} {100*(rl-PREM_RHO)/PREM_RHO:+11.1f}% "
              f"{100*(ro-PREM_RHO)/PREM_RHO:+11.1f}% {100*(ro-rl)/rl:+10.2f}%")
    print("\nThe own-derived ledger shifts the aggregate density by the 'own vs lit' column; light "
          "elements lower it in both. Wrote eos_params_qha.json.")


if __name__ == "__main__":
    main()
