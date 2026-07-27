#!/bin/bash
#SBATCH --job-name=qha_core
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --nodes=1                 # one node per job -- many small jobs, NOT one big job
#SBATCH --ntasks-per-node=32      # MPI ranks = cores; raise to 48 for the >96-atom cells
#SBATCH --cpus-per-task=1         # OMP=1 path
#SBATCH --time=06:00:00           # 96-160 atom phonon supercells: 2 h was too short (whole arrays TIMEOUT), 6 h has margin
#SBATCH --partition=public        # ASU public partition ('general' went private-only May 2026; NOT highmem)
# htc (4 h limit) is also good for these short jobs; do NOT request a highmem node
#SBATCH --mem=0                   # all RAM on the (standard) node -- supercell phonon SCFs need >16G
#
# Submit with the array + manifest chosen at submit time (see RUNBOOK.md), e.g.
#   STATIC:  sbatch --array=$(python idx.py static) --export=ALL,MANIFEST=job_manifest.csv  templates/slurm_array.sh
#   PHONON:  sbatch --array=0-$(($(wc -l < phonon_jobs.csv)-2))%50 --export=ALL,MANIFEST=phonon_jobs.csv templates/slurm_array.sh
: "${MANIFEST:=job_manifest.csv}"

# ---- the Cu-Zr lessons, baked in ----
export OMP_NUM_THREADS=1          # the fix that unstuck the Cu-Zr SCFs
export MKL_NUM_THREADS=1
ulimit -s unlimited
mkdir -p logs

# Quantum ESPRESSO from the qha conda env (pw.x), not a module -- no VASP license needed.
module load mamba/latest
source activate qha
# On Sol, `source activate` sets the prompt but NOT the PATH -- force the env's bin to the front,
# else `python` is the system python and ase/phonopy are missing.
export PATH="$HOME/.conda/envs/qha/bin:$PATH"
export QE_CMD="mpirun -np ${SLURM_NTASKS:-32} pw.x"        # conda QE: mpirun, not srun
export ESPRESSO_PSEUDO="$HOME/pseudo/earthcore"            # TODO: your .UPF dir (Fe/O/Si/S/C/Ni)

# One array index -> one pw.x run (static E-V or phonon displacement), per $MANIFEST.
"$HOME/.conda/envs/qha/bin/python" run_one.py --index "${SLURM_ARRAY_TASK_ID}" --manifest "${MANIFEST}"
