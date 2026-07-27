# QHA core endmembers — campaign README

First-principles quasiharmonic F(V,T) for the FeSiOCS inner-core paper's solid phases, sized
to the **ASU Sol cores-only / 5-8 h wall / no-deep-memory** budget. Read `PLAN.md` first.

**DFT engine: Quantum ESPRESSO (`pw.x`)**, run from the conda env (no VASP license). See `RUNBOOK.md`.

## Layout
```
PLAN.md            campaign design, scope, validation gates (read first)
RUNBOOK.md         exact Sol run sequence (start here to execute)
make_manifest.py   builds the (phase x volume x displacement) job list + budget summary
job_manifest.csv   the job list (generated; one row per SCF)  <- regenerate, don't hand-edit
qe_common.py       shared QE settings (ECUTWFC/ECUTRHO/KSPACING/DEGAUSS/MAGMOM) + input writer
run_one.py         per-array-index driver: static vc-relax (E-V) or phonon scf, via pw.x
gen_displacements.py  phonopy --qe displacements from the relaxed cells -> phonon_jobs.csv
fit_qha.py         force sets -> F(V,T) -> EoS (skeleton)
idx.py             prints static-job indices for the SLURM array
structures/        starter POSCARs (read by ASE, converted to QE input)
templates/slurm_array.sh   SLURM job-array (OMP=1, one node/job, no highmem; activates qha env)
```

## Budget (from `python make_manifest.py`, ESTIMATES until pre-flight)
486 jobs, ~4,755 core-hours, every job <0.85 h wall and <=7.3 GB. Fits cores-only / short-wall.
Fe7C3 dominates (279 jobs) — it's the highest-uncertainty phase and the biggest cell; gate it.

## Workflow on Sol
1. **Reuse the Cu-Zr Zentropy/dfttk stack.** Same OMP=1 lesson, same phonopy pipeline.
2. **Pre-flight (PLAN.md):** hcp-Fe, one volume, end-to-end. Confirm SCF converges OMP=1 on a
   standard node, peak mem < ceiling, wall < budget, phonon DOS clean. Write the MEASURED
   cores/mem/wall back into `make_manifest.py::estimate()` and regenerate the manifest.
3. **Validate hcp-Fe** QHA rho(P,300K) vs Dewaele 2006 (the benchmark) BEFORE running FeO/Fe7C3.
4. **Scale out** per phase via the SLURM array; collect forces with phonopy; build F(V,T).
5. **Fit + emit** `eos_params_qha.json` (`own_derived_parameter` mgd3 block per phase, literature
   value kept as `validation_comparison`), then drop into `../core/eos_params.json`.

## Hard rules (Sol budget)
- One node per job, many small jobs. Never a single big DFPT/AIMD job.
- `OMP_NUM_THREADS=1`; MPI ranks = cores (the Cu-Zr fix).
- Standard/public QOS only. **Do not request a highmem node.** If a job needs one, shrink the
  supercell instead.
- Phase B (insertion free energies for D(T)) is GATED on Phase A memory headroom. Liquid Fe-O
  AIMD is OUT.

## TODO before production (expert judgment — not guessable)
- Sol partition/QOS/account names in `slurm_array.sh`; `ESPRESSO_PSEUDO` pointed at your `.UPF` set.
- `ECUTWFC`/`ECUTRHO`/`KSPACING`/`DEGAUSS`/`MAGMOM` frozen after the hcp-Fe pre-flight (`qe_common.py`).
- Supply `POSCAR_Fe3S` / `POSCAR_Fe7C3` (structures/README.md); confirm supercells/atom counts.
