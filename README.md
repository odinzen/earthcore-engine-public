# earthcore-engine

Reproducibility package for a thermodynamic accounting of the light elements in Earth's core, with a
first-principles equation-of-state validation of the solid endmembers. This repository holds code and
data only; the associated manuscripts are forthcoming.

## Layout

- **`core/`**: the inner-core composition model. A multiphase iron-alloy aggregate (hexagonal iron
  with B2 FeSi, metallic FeO, the carbide Fe7C3, the sulfide Fe3S, and nickel) is built from a sourced
  Mie-Grueneisen-Debye equation-of-state ledger (`eos_params.json`), evaluated with BurnMan, and closed
  against the Preliminary Reference Earth Model (PREM) by a chemical-potential balance across the
  inner-core boundary (`partition_params.json`). `core_sweep_solver.py` does the property matching;
  `full_space_rerun.py --production` runs the four-element double-constraint search, the sensitivity
  scenarios and the Monte Carlo, and writes `production/`. `production_checks.py`, `icb_circularity.py`,
  `icb_velocity.py` and `verify_qha_eos.py` are the checks; `reproduce.sh` runs everything in order.
  See `core/README.md` for the result and the file map.
- **`qha_campaign/`**: the first-principles quasiharmonic equation-of-state campaign for the solid
  endmembers. Density functional theory (Quantum ESPRESSO, the PBE functional, projector augmented wave
  pseudopotentials) finite-displacement phonons give the Helmholtz free energy F(V,T) per phase.
  `dump_qha.py` exports the per-phase free-energy tables (`qha_*.json`); `fit_eos_qha.py` fits the
  equation of state and emits `eos_params_qha.json`; `qha_benchmark_figure.py` produces the validation
  figure. The exact settings and job sequence are in `reproduce_qha.md`.

## Reproduce

```
conda env create -f environment.yml
conda activate earthcore-engine

cd core && bash reproduce.sh          # every result table and figure
python ../qha_campaign/fit_eos_qha.py # first-principles equation of state fit
```

The first-principles calculations were run on the Sol supercomputer at Arizona State University;
`qha_campaign/reproduce_qha.md` records the settings and the Sol job sequence.

## License

Code is licensed under the MIT License; data, result tables, and figures under CC BY 4.0.
See [`LICENSE`](LICENSE).
