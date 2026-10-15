#!/usr/bin/env bash
# Reproduce every result, table and figure of the Fe-Ni-Si-O-C-S core composition study.
set -e
python core_sweep_solver.py              # Stage 1 property matching
python liquid_carbon_badro.py            # liquid carbon component (Badro 2014 Fig. 1)
python full_space_rerun.py --production  # Stage 2 search, sensitivity, Monte Carlo -> production/
python production_checks.py              # D_C and Xie-oxygen checks -> production/checks.txt
python icb_circularity.py                # oxygen-requirement audit
python icb_velocity.py                   # Vp, Vs, Poisson of the best-fit solid
python verify_qha_eos.py                 # calculated vs imported solid EoS
python fe_si_o_c_figures.py              # figures 1-4
python fe_si_o_c_partition_figures.py    # figures 5-10
python make_fig11.py                     # figure 11
python qha_benchmark_figure.py           # figure 12
python graphical_abstract.py
echo "Done. Results in production/, figures in figures/."
