# earthcore-engine

Reproducibility package for a thermodynamic accounting of the light elements in Earth's core, with a
first-principles equation-of-state validation of the solid endmembers. This repository holds code and
data only; the associated manuscripts are forthcoming.

## Layout

- **`core/`** — the inner-core composition model. A multiphase iron-alloy aggregate (hexagonal iron
  with B2 FeSi, metallic FeO, the carbide Fe7C3, the sulfide Fe3S, and nickel) is built from a sourced
  Mie-Grueneisen-Debye equation-of-state ledger (`eos_params.json`), evaluated with BurnMan, and closed
  against the Preliminary Reference Earth Model (PREM) by a chemical-potential balance across the
  inner-core boundary (`partition_params.json`). Solvers: `core_sweep_solver.py` (property matching),
  `icb_partition_sweep.py` (the double-constraint search), `icb_partition_uq.py` (uncertainty),
  `icb_partition_sensitivity.py`, `icb_velocity.py`, and `icb_circularity.py` (the circularity audit).
  The `*.csv` files are the result tables; `reproduce.sh` regenerates them. `verify_qha_eos.py`
  cross-checks the result on the first-principles equation of state.
- **`qha_campaign/`** — the first-principles quasiharmonic equation-of-state campaign for the solid
  endmembers. Density functional theory (Quantum ESPRESSO, the PBE functional, projector augmented wave
  pseudopotentials) finite-displacement phonons give the Helmholtz free energy F(V,T) per phase.
  `dump_qha.py` exports the per-phase free-energy tables (`qha_*.json`); `fit_eos_qha.py` fits the
  equation of state and emits `eos_params_qha.json`; `qha_benchmark_figure.py` produces the validation
  figure. The exact settings and job sequence are in `reproduce_qha.md`.
- **`magma/`** — a metal-silicate partitioning model for the oxygen provenance of the core
  (`oxygen_solver.py`, `continuous_accretion.py`).

## Reproduce

```
conda env create -f environment.yml
conda activate earthcore-engine

cd core && bash reproduce.sh          # composition sweep and data tables
python ../qha_campaign/fit_eos_qha.py # first-principles equation of state fit + validation
python verify_qha_eos.py              # cross-check the inner-core result on the calculated EoS
```

The first-principles calculations were run on the Sol supercomputer at Arizona State University;
`qha_campaign/reproduce_qha.md` records the settings and the Sol job sequence.

## License

Apache-2.0. See [`LICENSE`](LICENSE).
