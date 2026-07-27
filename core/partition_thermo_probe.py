"""NEGATIVE RESULT (kept for the record). Computing D_i in-repo from integrated
free energies fails: mgd3 minerals carry no formation reference (F_0=0), so
differencing host vs host does not cancel the arbitrary EoS baseline and D(P) comes
out unphysical by many orders of magnitude, whereas Alfe 2002 D is ~P-insensitive.
A defensible in-repo D needs ab initio insertion energies or a high-P CALPHAD
database, neither of which exists here. Production uses literature-anchored D.
Run: python partition_thermo_probe.py"""
import numpy as np, fe_si_o_c_model as base, icb_partition_model as M
R = 8.314462; PP = M.PP
C1 = PP["liquid_outer_core"]["components_at_130GPa"]; C2 = PP["liquid_outer_core"]["components_at_329GPa"]
P1, P2 = 130.0, 329.0
HOST = {"O": ("feo", 1.0, 1), "Si": ("fesi", 1.0, 1), "S": ("fe3s", 3.0, 1), "C": ("fe7c3", 7.0/3.0, 3)}
D_LIT = {"O": 0.02, "S": 0.75, "Si": 0.80, "C": 6.3e-5}
def _v(el, p): f = (p-P1)/(P2-P1); return (C1[el]["V"]+f*(C2[el]["V"]-C1[el]["V"]))*1e-6
def _g(k, p, t): m = base._PHASES[k]; m.set_state(p*1e9, t); return m.molar_gibbs
def _s(el, p, t): k, nfe, d = HOST[el]; return (_g(k, p, t)-nfe*d*_g("hcp_fe", p, t))/d
def _l(el, p, t, pr=328.85):
    _, nfe, _ = HOST[el]; ps = np.linspace(pr, p, 40)
    return np.trapz([(_v(el, q)-nfe*_v("Fe", q)) for q in ps], ps*1e9)
def D(el, p, t, pa=328.85, ta=5800.0):
    F = lambda P, T: _s(el, P, T)-_l(el, P, T); C = -R*ta*np.log(D_LIT[el])-F(pa, ta)
    return np.exp(-(F(p, t)+C)/(R*t))
if __name__ == "__main__":
    T = 5800.0
    print("Attempted in-repo D(P) anchored at ICB. UNPHYSICAL:")
    for el in ("O", "Si", "S", "C"):
        print(f"  {el}: D(135)={D(el,135.75,T):.2e}  D(329 anchor)={D(el,328.85,T):.2e}  lit={D_LIT[el]:.2e}")
