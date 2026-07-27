# Fe Ni Si O C S inner core and outer core against PREM

A thermodynamic accounting of the light elements in Earth's core. A sourced
equation of state ledger builds a multiphase iron alloy aggregate (close packed
hexagonal iron with B2 FeSi, metallic FeO, the carbide Fe7C3, the sulfide Fe3S,
plus nickel), matches its density and bulk sound velocity to PREM, and closes the
composition with a chemical potential balance across the inner core boundary
(Badro 2014 double constraint). Engine: BurnMan, Mie Grueneisen Debye EoS.

## Result

Admissible outer core: oxygen rich, 6 to 19 mol percent oxygen (best near 18, about
5.9 wt percent), one to two mol percent sulfur, silicon near zero, nickel 5 atomic
percent. Oxygen is required, recovered from the partition step and shown structural
by icb_circularity.py (no oxygen-free liquid fits; the requirement survives a 50
percent degradation of the oxygen calibration). Inner core: iron and nickel with
about one mol percent sulfur and trace oxygen and carbon. Density jump near 5.6
percent brackets the observed 4.5 to 4.8. Aggregate Vp matches PREM; Vs overshoots
by 28 percent, the inner core shear anomaly (premelting, not composition).

## Files

Ledgers: eos_params.json (solid phases), partition_params.json (Fe3S, liquid arm,
partition coefficients, nickel, PREM outer core, shear moduli). Model:
fe_si_o_c_model.py (solid aggregate, PREM, Fe thermal correction),
icb_partition_model.py (quinary solid with nickel, liquid arm, partition step).
Solvers: core_sweep_solver.py (Stage 1), icb_partition_sweep.py (Stage 2; pass
'carbide' for the carbon-into-solid scenario). Analyses: icb_velocity.py (Vp/Vs/
Poisson), icb_partition_uq.py (Monte Carlo), icb_partition_sensitivity.py,
icb_circularity.py (referee circularity audit), partition_thermo_probe.py (documented
negative result on in-repo D). Figures: fe_si_o_c_figures.py, fe_si_o_c_partition_figures.py,
make_fig11.py, graphical_abstract.py. Docs: CIRCULARITY_AUDIT.md, the manuscript and
supplementary, environment.yml.

## Run (miniforge, Ubuntu)

    conda activate earthcore-engine
    python core_sweep_solver.py
    python icb_partition_sweep.py
    python icb_velocity.py
    python icb_circularity.py

Switches: FEO_PHASE=B1 selects the stable oxide; FE_THERMAL=0 disables the iron
thermal correction.

## Key references

Dziewonski & Anderson 1981 (PREM). Badro et al. 2014 PNAS (double constraint).
Alfe et al. 2002 EPSL (partition coefficients). Li et al. 2019 (carbon). McDonough
2003 (nickel). Dewaele 2006 (iron). Fischer 2014/2011 (FeSi, FeO). Mookherjee/
Nakajima 2011 (Fe7C3). Seagle 2006, Kamada 2014 (Fe3S). Kuwayama 2020, Ichikawa
2014 (liquid iron). Murphy 2013, Ono 2013, Martorell 2013 (shear). Myhill 2023 (BurnMan).
