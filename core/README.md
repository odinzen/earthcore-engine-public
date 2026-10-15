# Fe Ni Si O C S inner core and outer core against PREM

A thermodynamic accounting of the light elements in Earth's core. A sourced
equation of state ledger builds a multiphase iron alloy aggregate (close packed
hexagonal iron with B2 FeSi, metallic FeO, the carbide Fe7C3, the sulfide Fe3S,
plus nickel), matches its density and bulk sound velocity to PREM, and closes the
composition with a chemical potential balance across the inner core boundary
(the Badro 2014 double constraint). Engine: BurnMan, Mie Grueneisen Debye EoS.

## Result

With carbon, silicon and sulfur held to published limits (1.0, 2.4 and 1.9 wt%),
the outer core needs 2.2 to 3.6 wt% oxygen; the best fit is 2.6 wt% oxygen, 2.3 wt%
silicon, 1.6 wt% sulfur and 1.0 wt% carbon. Oxygen is required in every sensitivity
scenario and in 179 of 198 admissible Monte Carlo draws. The inner core cannot reach
PREM on oxygen, silicon and sulfur alone because oxygen stays in the liquid; it needs
carbon entering the solid (carbon partition coefficient 1.5) and holds about 1.4 wt%
carbon, 1.7 wt% silicon and 1.2 wt% sulfur. At equal partitioning or at the
first-principles coefficient of Li et al. 2019 no composition fits. The boundary
density jump is 7.2 percent at the best fit (5.8 to 7.2 percent across the Monte Carlo).

## Files

Ledgers: eos_params.json (solid phases), partition_params.json (Fe3S, liquid
components, partition coefficients, nickel, PREM outer core, shear moduli).
Liquid components: oxygen and sulfur solved from the Fe-S-O liquids of Badro 2014
Table S1; carbon from Badro 2014 Fig. 1 (digitize_badro_fig1.py reads the figure from
the published PDF, liquid_carbon_badro.py derives the component into
liquid_carbon_badro.json); silicon from Jin 2026 (liquid_refit_independent_data.py,
which also holds the Xie 2025 oxygen check).
Model: fe_si_o_c_model.py (solid aggregate, PREM, Fe thermal correction),
icb_partition_model.py (solid with nickel, liquid arm, partition step).
Stage 1: core_sweep_solver.py. Stage 2: full_space_rerun.py --production runs the
four element search, the eight sensitivity scenarios and the 400 draw Monte Carlo,
and writes production/ (sweeps, sensitivity_summary.csv, uq_draws.csv,
uq_summary.json, best_fit.json). The Monte Carlo replays the draws recorded in
icb_partition_uq_results.csv by the original icb_partition_uq.py, and checks the
replay before it runs.
Checks: production_checks.py (carbon coefficient 1.0, 2.0 and Li 2019; Xie oxygen),
icb_circularity.py (oxygen-free search and degraded oxygen), icb_velocity.py
(Vp, Vs, Poisson), verify_qha_eos.py (calculated versus imported solid EoS),
partition_thermo_probe.py (documented negative result on calculating D here).
Figures: fe_si_o_c_figures.py (1-4), fe_si_o_c_partition_figures.py (5-10),
make_fig11.py, qha_benchmark_figure.py (12), graphical_abstract.py; all read
production/best_fit.json for the best fit composition.

## Run

    conda activate earthcore-engine
    bash reproduce.sh

Switches: FEO_PHASE=B1 selects the stable oxide; FE_THERMAL=0 disables the iron
thermal correction.

## Key references

Dziewonski & Anderson 1981 (PREM). Badro et al. 2014 PNAS (double constraint,
liquid components). Kuwayama et al. 2020 PRL (liquid iron). Jin et al. 2026
Minerals (liquid Fe-Si). Xie et al. 2025 (liquid Fe-O). Alfe et al. 2002 EPSL
(partition coefficients). Li et al. 2019; Chen et al. 2014; Prescher et al. 2015
(carbon). McDonough 2003 (nickel, sulfur and carbon limits). Dewaele 2006 (iron).
Fischer 2014/2011 (FeSi, FeO). Mookherjee 2011 (Fe7C3). Seagle 2006, Kamada 2014
(Fe3S). Ichikawa 2014 (liquid iron). Murphy 2013, Ono 2013, Martorell 2013 (shear).
Myhill 2023 (BurnMan).
