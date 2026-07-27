"""Print the comma-separated manifest indices for a given job_kind, for SLURM --array.

    sbatch --array=$(python idx.py static) ... templates/slurm_array.sh

(static_EV jobs are interleaved in job_manifest.csv, so an explicit index list is needed.)
"""
import csv
import os
import sys

kind = (sys.argv[1] if len(sys.argv) > 1 else "static").replace("static", "static_EV")
man = os.path.join(os.path.dirname(os.path.abspath(__file__)), "job_manifest.csv")
idx = [str(i) for i, r in enumerate(csv.DictReader(open(man))) if r["job_kind"] == kind]
print(",".join(idx))
