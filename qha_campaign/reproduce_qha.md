# Reproducing the QHA core-endmember EoS campaign

First-principles quasiharmonic F(V,T) for the FeSiOCS inner-core solid endmembers, feeding the
own-derived equation of state in `eos_params_qha.json`. This is the method record for the paper's
Methods section and for reproduction. Runs on ASU Sol (cores-only, no VASP license).

## Codes and pseudopotentials
- Quantum ESPRESSO **7.5** (`pw.x`), phonopy **4.2.1**, ASE for structure I/O.
- Pseudopotentials: **PBE (GGA), PAW, pslibrary 1.0.0** (Dal Corso, Comput. Mater. Sci. 95, 337, 2014).
  Fe and Ni use the semicore spin-polarized variants `Fe.pbe-spn-kjpaw_psl.1.0.0.UPF`,
  `Ni.pbe-spn-kjpaw_psl.1.0.0.UPF`; C, O, S, Si use `X.pbe-n-kjpaw_psl.1.0.0.UPF`.

## DFT settings (frozen after the hcp-Fe pre-flight; qe_common.py)
- ecutwfc 60 Ry, ecutrho 480 Ry (8x).
- Monkhorst-Pack k-mesh at KSPACING 0.20 Ang^-1 (2pi convention); the supercell mesh is the unit-cell
  mesh folded by the supercell dimension. Irreducible k-points per phase (static primitive cell):
  hcp-Fe 546, B2-FeSi 120, B1-FeO 145, Fe3S 100, Fe7C3 40.
- Methfessel-Paxton smearing, degauss 0.02 Ry.
- Spin-polarized (nspin=2) for the Fe-O/S/C phases; nspin=1 for hcp-Fe and B2-FeSi (nonmagnetic at
  core compression, and nspin=2 supercells magnetically diverge there). Starting moment Fe 2.2, Ni 0.6.
- mixing_mode = local-TF, mixing_beta 0.2, electron_maxstep 300 (stability for the magnetic supercells).
- conv_thr 1e-8 (static), 1e-7 (phonon displacements).

## Volume and temperature grids
- 9 phonon volumes at V/V0 = 0.55-0.78 (brackets ~100-400 GPa: CMB 136 to center 364) carrying F(V,T);
  3 static-only anchor volumes at V/V0 = 0.86, 0.94, 1.02 added to reach P~0 and pin ambient V0/K0/K0'.
- Static E-V: vc-relax at fixed volume (cell_dofree = 'shape'; shape + ions relax, V held at V0*ratio).
- Phonons: phonopy finite displacements, distance 0.01 Ang; supercells hcp 4x4x3 (96 at), FeSi 4x4x4
  (128), B1-FeO 3x3x3 (54), Fe3S 2x2x2 (128), Fe7C3 2x2x2 (160). Forces -> force constants ->
  phonon DOS -> F_vib(V,T) on a 300-7000 K grid (T is post-processing, no extra DFT).

## Sol sequence
1. `python make_manifest.py`  -> job_manifest.csv (phase x volume x displacement).
2. Static array: `sbatch --array=$(python idx.py static) --export=ALL,MANIFEST=job_manifest.csv templates/slurm_array.sh`.
3. `python gen_displacements.py`  -> phonopy displacements + phonon_jobs.csv.
4. Phonon array: `sbatch --array=0-$(($(wc -l < phonon_jobs.csv)-2))%50 --export=ALL,MANIFEST=phonon_jobs.csv templates/slurm_array.sh` (12 h walltime; run_one.py skips finished jobs via .done markers).
5. `python dump_qha.py --phase <P>`  -> qha_<P>.json (V, E_static, F_vib(T)); download via OnDemand.

## Off-cluster fit (numpy/scipy only; verifiable, no phonopy)
`python fit_eos_qha.py`:
- Cold isotherm: Birch-Murnaghan 3rd order to the static E(V) over V/V0 = 0.55-1.02 -> V0, K0, K0'.
- Thermal: direct thermodynamic Grueneisen gamma = V (dS/dV)_T / C_V from F(V,T), high-T average,
  outliers trimmed -> gamma0, q; effective Debye theta0 with gamma0/q fixed.
- Density gate: fit-free -dF/dV at 300 K vs the imported literature EoS across the core range.
- Emits eos_params_qha.json (own mgd3 per phase + literature validation_comparison + reliability flag).

## Phase notes
- **hcp-Fe** is the benchmark; density reproduces Dewaele 2006 to <=1.3% across the core. Its ambient
  K0 (nonmagnetic PBE, ~313 GPa) differs from Dewaele's 163 by construction (their fit runs through
  the magnetic collapse); the two agree in the core where both are nonmagnetic.
- **B1-FeO is dynamically unstable at core compression** (imaginary phonon modes, -7 to -43 THz, at
  every volume; not a mesh/convergence artifact). Its cold EoS (from the static energies) is own and
  validated, but QHA cannot give its vibrational thermodynamics, so the FeO Grueneisen term uses the
  literature value. The metallic B8-FeO polymorph would be the phase to compute for an own FeO thermal.
- **Fe3S, Fe7C3** own-derived cold + thermal; Fe7C3 (the highest-uncertainty imported phase) agrees
  with the literature in both density (<=1.9%) and gamma0 (1.62 vs 1.70).

Cross-check that the oxygen-rich core survives on the own EoS: `../core/verify_qha_eos.py` (needs burnman).
