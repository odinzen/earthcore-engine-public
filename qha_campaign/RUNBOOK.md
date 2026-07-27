# RUNBOOK — running the QHA campaign on ASU Sol

Step-by-step to take this from a plan to own first-principles F(V,T) for the core endmembers.
Read `PLAN.md` for the why. **You drive the Sol runs (login/2FA); this is the exact sequence.**
Hard rules (the Sol budget): cores-only, 5-8 h wall/job, **no highmem node**, `OMP_NUM_THREADS=1`.

**DFT engine: Quantum ESPRESSO (`pw.x`)** — installed into the conda env, no VASP license / module
needed. phonopy drives it natively.

## 0. One-time setup on Sol
```bash
git clone https://github.com/odinzen/earthcore-engine.git    # (token auth; see session notes)
cd earthcore-engine/qha_campaign

# env with phonopy + ase + Quantum ESPRESSO (all conda-forge, no license)
module load mamba/latest
mamba create -n qha -c conda-forge python=3.10 phonopy ase qe numpy -y
source activate qha
which pw.x                                  # confirm QE installed

# free pseudopotentials (SSSP PBE efficiency) -> a dir of .UPF files
mkdir -p ~/pseudo/sssp_pbe_efficiency && cd ~/pseudo/sssp_pbe_efficiency
#   download SSSP PBE efficiency from materialscloud.org/discover/sssp (tar of .UPF),
#   or any QE .UPF set (PseudoDojo / PSLibrary). Unpack so Fe.*.UPF, O.*.UPF, ... sit here.
cd -
export QE_CMD="srun pw.x"
export ESPRESSO_PSEUDO="$HOME/pseudo/sssp_pbe_efficiency"
```
Edit `templates/slurm_array.sh`: set your **partition / QOS / account** (standard/public — NOT
highmem) and the `ESPRESSO_PSEUDO` path, then `python make_manifest.py` to (re)write the manifest.

## 1. Supply the two complex structures
`structures/` ships starter POSCARs for hcp_Fe, B2_FeSi, B1_FeO, B8_FeO. **You add**
`POSCAR_Fe3S` and `POSCAR_Fe7C3` (see `structures/README.md`), confirm their atom counts match
`make_manifest.py`, regenerate the manifest if you changed any cell.

## 2. PRE-FLIGHT — hcp-Fe, one volume, end to end (DO THIS BEFORE ANYTHING SCALES)
Pick a mid-compression hcp_Fe static job and run it interactively:
```bash
# index of hcp_Fe static at v04 (mid grid)
python -c "import csv;print([i for i,r in enumerate(csv.DictReader(open('job_manifest.csv'))) if r['phase']=='hcp_Fe' and r['job_kind']=='static_EV' and r['vol_idx']=='4'][0])"
salloc --nodes=1 --ntasks=32 --time=1:00:00 --partition=public   # cores-only, no highmem ('general' is private since May 2026)
source activate qha
export OMP_NUM_THREADS=1 QE_CMD="srun pw.x" ESPRESSO_PSEUDO="$HOME/pseudo/sssp_pbe_efficiency"
python run_one.py --index <that-index>
```
Confirm, and only then proceed: (a) SCF converged with OMP=1 on a standard node, (b) peak memory
< node ceiling, (c) wall < budget, (d) `ECUTWFC`/`ECUTRHO`/`KSPACING`/`DEGAUSS` are converged (bump
each, watch total energy / pressure settle), (e) magnetism is sane per phase.
**Freeze those values in `qe_common.py`** (they're shared by the static and phonon inputs). Record
the measured cores/mem/wall back into `make_manifest.py::estimate()`.

**Validation gate:** check hcp-Fe rho(P, 300 K) from the relaxed `static/pw.out` against Dewaele
2006 (V0=6.753 cm3/mol-atom, K0=163.4 GPa, K'=5.38) to within a few percent. If it misses, the
settings (pseudopotential, ECUTWFC, k-grid, smearing) are wrong — fix before trusting FeO/Fe7C3.

## 3. Stage 1 — static E(V) array (54 jobs: 6 phases x 9 volumes)
```bash
sbatch --array=$(python idx.py static)%50 --export=ALL,MANIFEST=job_manifest.csv templates/slurm_array.sh
```
Each writes `runs/<phase>/v<idx>/static/pw.out` (the relaxed cell + E(V) at that volume).

## 4. Stage 2 — generate displacements (cheap; login node or a tiny job)
```bash
python gen_displacements.py          # reads every static pw.out; writes disp/supercell-*.in + phonon_jobs.csv
```
This writes `phonon_jobs.csv` — the **real** phonon job list (phonopy's actual displacement count,
which may differ from the manifest estimate).

## 5. Stage 3 — phonon force array (driven off the real list)
```bash
N=$(($(wc -l < phonon_jobs.csv) - 1))
sbatch --array=0-$((N-1))%50 --export=ALL,MANIFEST=phonon_jobs.csv templates/slurm_array.sh
```
Gate order: let **hcp-Fe finish and validate first**, then FeO (oxygen arm), and **Fe7C3 last** —
highest uncertainty + biggest cell (most jobs). Use `--phase` filtering by editing the array if you
want to stage phase-by-phase.

## 6. Post-process — own F(V,T) -> EoS
```bash
python fit_qha.py --phase hcp_Fe     # then FeO, FeSi, Fe3S, Fe7C3
```
Builds FORCE_SETS, thermal properties, then (TODO in the skeleton) fits F(V,T) -> mgd3 params ->
`eos_params_qha.json`, keeping the imported literature value as `validation_comparison`.

## 7. Feed back into the paper
Drop each phase's own block into `../core/eos_params.json` tagged `own_derived_parameter` (keep the
sourced value beside it). Re-run `../core/reproduce.sh`: the headline numbers should move within the
stated uncertainty, and the manuscript's §2.1 / §5 "companion first-principles campaign" stops being
a promise and becomes a result.

## Troubleshooting (Cu-Zr lessons)
- SCF won't converge / OOM -> too-big/too-memory job. Keep `OMP=1`, shrink the supercell, or lower
  `KSPACING` density / `ECUTWFC`. Never request a highmem node — that's the delayed path.
- `no .UPF for element X` -> `ESPRESSO_PSEUDO` doesn't hold a pseudopotential for X; add it.
- A whole phase's static jobs fail identically -> bad starter POSCAR or a missing pseudopotential;
  check `structures/` and `$ESPRESSO_PSEUDO`.
- Phonon DOS has large imaginary branches -> cell not actually relaxed at that volume, or supercell
  too small; re-relax the static cell, enlarge `--dim`.
