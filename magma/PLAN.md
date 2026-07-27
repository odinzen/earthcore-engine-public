# Oxygen provenance: where did the core's oxygen come from?

## Question
Our core model needs ~5 wt% O in the outer core. Find the deep magma-ocean state
(depth/P, T, oxygen fugacity / mantle FeO) that simultaneously:
  (a) drives ~5 wt% O into the metal via sourced metal-silicate O partitioning,
  (b) reproduces mantle siderophile abundances (Ni, Co, Cr, V),
  (c) matches PREM outer-core density.
Then: the leftover mantle FeO predicts the chemistry of the magmas it makes.

## Novelty vs prior art
Badro 2015 (PNAS) already did the coupled inversion (2.7-5 wt% O). Our delta:
start from a FIXED ~5 wt% O (high end, vs Rubie ~0.5), carry full uncertainty
quantification, and add the falsifiable mantle-magma prediction.

## Method (all literature-anchored, no new DFT/MD)
- O partitioning D_O(T, P, fO2): Rubie 2004, Ricolleau 2011, Fischer 2015.
- Siderophile match (Ni,Co,Cr,V): Wood 2006, Siebert 2013, Fischer 2015.
- Self-consistency solver + Monte Carlo: reuse the core pipeline's UQ machinery.
- Mantle-magma prediction: Katz 2003 / pyMelt (optional, lightweight).

## Cost vs impact
Compute: CHEAP. Algebra + Monte Carlo, laptop seconds to minutes. No DFT/MD.
Ingredients (D values, EoS) are published; we ingest, not compute.
Impact: solid follow-on. Realistic target PEPI / EPSL / GCA (not PNAS as a sequel).
Verdict: very favorable cost/impact ratio; high synergy with the core paper.

## Deliverables
ledger (sourced D_O + siderophile params), solver, UQ, 3-4 figures, manuscript.
