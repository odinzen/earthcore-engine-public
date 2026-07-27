# QHA core endmembers — first-principles F(V,T) for the FeSiOCS inner-core paper

## Why this exists
The FeSiOCS core paper (`../core/`) currently builds every solid phase from
**imported** Mie-Grueneisen-Debye parameters (`eos_params.json`, Dewaele/Fischer/Mookherjee).
That is honest but it is *accounting of other people's EoS numbers* — a scope risk for the
retarget to **Matter & Radiation at Extremes (MRE)**, which wants first-principles
extreme-conditions physics.

This campaign computes **our own quasiharmonic (QHA) Helmholtz free energy F(V,T)** for the
solid endmembers at core P-T, so the paper reads as *first-principles high-pressure
thermodynamics of the core* rather than a literature ledger. It is the **load-bearing upgrade
content** that earns the venue — not an optional add-on.

## What it produces (and where it lands in Paper 1)
1. **hcp-Fe** — QHA F(V,T) validated against Dewaele 2006 (the benchmark, mirroring the
   Cu-Zr B0-vs-Du-2014 check). This is the credibility anchor.
2. **B8/B1-FeO** — own F(V,T) for the **oxygen arm**. The headline result is the oxygen
   requirement; replacing the imported FeO thermal terms with own first-principles values
   moves the "calibrated proxy" critique off the *solid* side.
3. **Fe7C3** — own F(V,T) for the **highest-uncertainty phase** (`eos_params.json` FLAG:
   thermal MGD calibrated to ~30 GPa and extrapolated ~10x to ICB; shear calibrated to
   velocity, not independent). Retiring this is the single biggest scientific win.
4. **B2-FeSi, Fe3S** — own F(V,T) to complete the ledger so the whole solid aggregate is
   first-principles, not mixed-provenance.

Each F(V,T) is fit back to the `mgd3` form and **replaces** the corresponding block in
`eos_params.json`, tagged `own_derived_parameter` (literature-only / own-derived discipline).
Imported values are kept alongside as the validation comparison, not deleted.

## Budget constraints (ASU Sol — HARD)
- **Cores only.** No deep-memory / highmem nodes (the Cu-Zr SCF stalls were memory + convergence).
- **5-8 h wall per job, max.** Structure the campaign as **many small independent jobs**, not
  a few large ones.
- **OMP_NUM_THREADS=1** (the Cu-Zr fix). Parallelize over k-points/bands with MPI ranks = cores.
- Standard/public partition QOS; SLURM **job array** so hundreds of jobs drain through the
  per-user limit without babysitting.

These constraints *fit* QHA naturally: finite-displacement (supercell) phonons decompose into
many independent static SCFs — the opposite of one big DFPT/AIMD job. We deliberately use
**finite-displacement phonons over DFPT** because each displacement is a separate short,
standard-memory SCF.

## Method
For each (phase, volume):
1. Static SCF relax of the primitive cell at fixed V (ISIF=2-style; volume set, ions+shape relaxed).
2. `phonopy` finite-displacement supercell → N_disp independent SCF force jobs.
3. Force sets → phonon DOS → QHA F(V,T) on a post-hoc T-grid (cheap, no extra DFT).
Then per phase: fit F(V,T) across the volume grid → Vinet/BM3 (cold) + Debye-Grueneisen (thermal)
→ emit `mgd3` parameters.

- **Volume grid:** 9 points spanning V/V0 in [0.55, 0.78] per phase (brackets ~100-400 GPa,
  i.e. CMB 136 to center 364 GPa, with margin to fit the EoS).
- **T-grid (post-processing only):** 300-7000 K, covering the 5500-6200 K ICB scan plus margin.
- **DFT jobs = phases x volumes x N_disp** only (T is free). See `make_manifest.py` for the count.

### Phase A vs Phase B
- **Phase A — endmember QHA (this campaign, fits the budget).** The 5-6 phases above. Small
  cells, standard memory, short wall. Fe7C3 and Fe3S have the largest supercells and are the
  ones that approach the ceiling — use the smallest symmetry-sufficient supercell and split the
  displacement array; if a single Fe7C3 displacement SCF exceeds the memory/wall ceiling in the
  pre-flight, drop the supercell size before scaling out.
- **Phase B — insertion free energies for D(T) (GATED).** The temperature-dependent partition
  coefficients that the in-repo derivation failed to get (`partition_thermo_probe.py`, manuscript
  §2.6) need *defected* supercells (one light atom in hcp-Fe) — bigger cells, more memory. Only
  attempt after Phase A confirms the per-job memory headroom. Do **not** let Phase B block Phase A.
- **OUT OF SCOPE:** liquid Fe-O AIMD (the Badro-independent oxygen test). Long wall + deep
  memory. Stays future work, per the Sol budget.

## Pre-flight (do this before scaling out — 1 phase, 1 volume)
Run hcp-Fe at one core-compression volume end-to-end: confirm (a) SCF converges with OMP=1 on a
standard node, (b) per-job peak memory < node ceiling, (c) wall < budget, (d) the phonon DOS is
gap-free / no large imaginary branch. Record the measured cores/mem/wall back into
`make_manifest.py` so the full-campaign budget estimate is real, not assumed.

## Validation gates (no result ships without these)
- hcp-Fe QHA rho(P, 300 K) within a stated % of Dewaele 2006 Vinet (target: a few %, like the
  Cu-Zr B0 check). If it misses, the settings are wrong — fix before trusting FeO/Fe7C3.
- Each phase's QHA-fit K_0, K'_0 compared to the imported literature value in `eos_params.json`;
  agreement = confidence, disagreement = a documented finding (and a manuscript talking point).

## Deliverables
- `eos_params_qha.json` — own first-principles `mgd3` block per phase + the literature value kept
  as `validation_comparison`.
- A validation table (QHA vs imported) → a manuscript table + the hcp-Fe benchmark figure.
- `reproduce_qha.md` — exact Sol commands, settings, SLURM array script, phonopy version.
- Two manuscript paragraphs: the first-principles method, and the validation/benchmark.

## Stack
Quantum ESPRESSO (`pw.x`) + `phonopy`, run from the conda env (no VASP license; VASP is gated on
Sol). `qe_common.py` holds the shared QE settings; see `make_manifest.py` for the job matrix and
`templates/slurm_array.sh` for the array (cores-only, OMP=1, no-highmem QOS). Full sequence: `RUNBOOK.md`.
