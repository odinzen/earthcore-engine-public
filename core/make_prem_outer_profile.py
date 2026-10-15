"""Tabulate PREM's outer core (Dziewonski and Anderson 1981) against pressure for Fig. 5.

The fit uses PREM only at its two boundaries (prem_outer at the CMB and ICB); this table
lets the figure draw the real profile in between instead of a straight line.
Run in the project environment (BurnMan 2.1.0).
"""
import numpy as np
import burnman

r_km = np.linspace(1221.5 + 0.5, 3480.0 - 0.5, 120)
P, rho, vp = burnman.seismic.PREM().evaluate(["pressure", "density", "v_p"], (6371.0 - r_km) * 1000.0)
order = np.argsort(P)
np.savetxt("prem_outer_core_profile.csv",
           np.c_[P[order] / 1e9, rho[order] / 1000.0, vp[order] / 1000.0], delimiter=",", fmt="%.5f",
           header="PREM outer core (Dziewonski and Anderson 1981) via BurnMan 2.1.0 seismic.PREM; "
                  "pressure_GPa,density_g_cm3,vphi_km_s (vs=0 so vphi=vp)", comments="# ")
