# Referee critique: audit against the code base

A simulated PEPI review raised several objections. Each is audited here against the
actual pipeline, and the decisive one (circularity of the oxygen requirement) is
settled with a direct computation (icb_circularity.py).

## Claim by claim

1. Liquid model is calibrated to Badro 2014. TRUE. Component volumes and moduli in
   partition_params.json are solved to reproduce Badro's ab initio alloy points.
   Audited and addressed by the circularity test below.

2. Partition coefficients are imported, not derived here. TRUE and stated. The
   in-repo derivation fails (partition_thermo_probe.py): mgd3 carries no formation
   reference, so D(P) comes out unphysical. The paper frames the work as synthesis.

3. Inner core matched on density and Vphi, not Vs. TRUE. icb_velocity.py: the
   best-fit solid gives Vp 11.45 (+3.8 percent), Vphi 10.21 (match), but Vs 4.48
   versus PREM 3.50 (+27.9 percent) and Poisson 0.41 versus 0.44. This is the inner
   core shear anomaly; the resolution is premelting and anelastic softening, not
   composition. Now reported as a result, section 3.5.

4. Stack of approximations (F_0=0, Fe7C3 extrapolated, FeO shear gap, single T).
   TRUE and each is flagged in the ledger; the uncertainty quantification and the
   sensitivity analysis bound their effect.

## The decisive test: is the oxygen requirement circular?

Run: python icb_circularity.py

Test A. Can ANY oxygen-free liquid satisfy the double constraint? Result: zero
admissible oxygen-free models. The best oxygen-free liquid (Si 12, S 0.5, C 4 mol
percent) FAILS on outer-core DENSITY at both the core mantle boundary and the inner
core boundary. Silicon, sulfur, and carbon alone cannot close the outer-core density
deficit. The oxygen requirement is therefore structural, set by PREM, not by the
oxygen calibration.

Test B. Degrade the oxygen component off its Badro fit (move its partial molar
volume and modulus toward iron, weakening its density effect) and re-run.

  O potency 1.0 (Badro):  admissible 45,  oxygen 10 to 19 mol percent,  required YES
  O potency 0.7:          admissible 11,  oxygen 9 to 14 mol percent,   required YES
  O potency 0.5:          admissible 3,   oxygen 9 to 12 mol percent,   required YES
  O potency 0.3:          admissible 0,   no solution,                  required (n/a)
  O potency 0.0:          admissible 0,   no solution,                  required (n/a)

Oxygen survives a 50 percent degradation of its calibration and stays required. As
its potency falls further the admissible window shrinks and then vanishes, but at no
point does an admissible low-oxygen or oxygen-free model appear. The requirement is
not an artifact of the precise calibration.

## Verdict

The circularity charge does not hold: the oxygen requirement is driven by the PREM
density constraint, which only oxygen among the candidates can satisfy at its
sourced effectiveness, and it is robust to substantial degradation of the oxygen
calibration. The imported partition coefficients and the calibrated liquid proxy
remain honest limitations, stated plainly. The inner core is matched in density and
bulk sound velocity; the shear velocity overshoot is the known anomaly, reported and
attributed to anelastic softening rather than composition.
