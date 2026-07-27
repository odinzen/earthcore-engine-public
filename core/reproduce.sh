#!/usr/bin/env bash
# Reproduce all results and data tables for the Fe-Ni-Si-O-C-S core composition study.
set -e
python core_sweep_solver.py
python icb_partition_sweep.py
python icb_partition_sweep.py carbide
python icb_velocity.py
python icb_circularity.py
echo "Done. See *.csv and CIRCULARITY_AUDIT.md."
