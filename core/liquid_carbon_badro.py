"""Liquid carbon component from Badro et al. (2014) Fig. 1, in place of "carbon equal to oxygen".

Badro et al. computed liquid Fe(1-x)C(x) and Fe(1-x)O(x) at x = 8.3 and 16.7 mol% at CMB and ICB
conditions with the same ab initio setup, but published those binaries only as plotted points
(Fig. 1); their Table S1 holds pure Fe and two Fe-S-O compositions, which is what the production O
and S components reproduce. The points below were digitized from the Fig. 1 raster (PNAS 111, 7542,
page 2) by colour-blob centroid after calibrating each axis on its tick marks; one pixel is about
5 kg/m3 in density and 2.3 m/s in velocity. The calibration reproduces Badro's tabulated pure Fe at
the ICB (12897 vs 12901 kg/m3, 10195 vs 10194 m/s). Where C and O markers overlap (ICB density, as
Badro's caption states), both take the visible marker. The CMB velocity C marker at 8.3 mol% sits
under the dashed target line and is read at the line (occluded, about +/-10 m/s).

Carbon enters the proxy as the ratio of Badro's Fe-C to Fe-O at 16.7 mol%, applied to the
proxy's own Fe-O, so the proxy's pure Fe anchors and Badro's shared systematics both cancel. The
8.3 mol% points are held out as a check on ideal mixing.

Run:  python liquid_carbon_badro.py   (writes liquid_carbon_badro.json)
"""
import json
import os

import icb_partition_model as M

HERE = os.path.dirname(os.path.abspath(__file__))

# Badro 2014 Fig. 1, digitized. (rho kg/m3, vphi m/s) keyed by boundary and mol% light element.
FIG1 = {
    "CMB": {"Fe": (10617, 7949),
            "O": {8.3: (10297, 8035), 16.7: (9959, 8100)},
            "C": {8.3: (10338, 8058), 16.7: (10039, 8151)}},
    "ICB": {"Fe": (12897, 10195),
            "O": {8.3: (12554, 10275), 16.7: (12177, 10341)},
            "C": {8.3: (12554, 10348), 16.7: (12177, 10444)}},
}
ANCHOR_P = {"CMB": 130.0, "ICB": 329.0}   # the proxy's two calibration pressures
X_FIT, X_CHECK = 16.7, 8.3


def binary(el, x, p):
    """Proxy density and bulk sound velocity of Fe(1-x)X(x), no nickel."""
    kw = dict(x_c=0.0, x_o=0.0, x_si=0.0, x_s=0.0)
    kw["x_c" if el == "C" else "x_o"] = x
    return M.liquid_props(kw["x_c"], kw["x_o"], kw["x_si"], kw["x_s"], p, x_ni=0.0)


def solve_carbon(boundary):
    p = ANCHOR_P[boundary]
    x = X_FIT / 100
    rho_o, v_o = binary("O", x, p)
    (rb_o, vb_o), (rb_c, vb_c) = FIG1[boundary]["O"][X_FIT], FIG1[boundary]["C"][X_FIT]
    rho_c, v_c = rho_o * rb_c / rb_o, v_o * vb_c / vb_o           # target Fe-C in proxy units
    v_fe, k_fe = M._comp("Fe", p)
    m_mix = (1 - x) * M.M_ATOM["Fe"] + x * M.M_ATOM["C"]
    v_mix = m_mix / (rho_c * 1000) * 1e6                            # cm3/mol
    k_mix = (v_c * 1000) ** 2 * rho_c * 1000 / 1e9                  # GPa
    V_C = (v_mix - (1 - x) * v_fe) / x
    K_C = x / (1 / k_mix - (1 - x) / k_fe)
    return V_C, K_C


def main():
    comps = {b: solve_carbon(b) for b in ("CMB", "ICB")}
    liq = M.PP["liquid_outer_core"]
    old = {b: dict(liq[f"components_at_{int(ANCHOR_P[b])}GPa"]["C"]) for b in comps}
    for b, (V, K) in comps.items():
        liq[f"components_at_{int(ANCHOR_P[b])}GPa"]["C"].update(V=V, K=K)

    print("carbon component (V cm3/mol, K GPa):")
    for b, (V, K) in comps.items():
        o = old[b]
        print(f"  {b} {ANCHOR_P[b]:.0f} GPa: V {V:.3f} (was {o['V']:.3f}), K {K:.1f} (was {o['K']:.1f})")
    print("held-out check at 8.3 mol% (proxy C/O ratio vs Badro C/O ratio):")
    for b in ("CMB", "ICB"):
        p, x = ANCHOR_P[b], X_CHECK / 100
        ro, vo = binary("O", x, p); rc, vc = binary("C", x, p)
        (rbo, vbo), (rbc, vbc) = FIG1[b]["O"][X_CHECK], FIG1[b]["C"][X_CHECK]
        print(f"  {b}: rho {100*(rc/ro-1):+.2f}% vs {100*(rbc/rbo-1):+.2f}%,"
              f"  vphi {100*(vc/vo-1):+.2f}% vs {100*(vbc/vbo-1):+.2f}%")

    out = {"_README": __doc__.split("\n\n")[0],
           "components": {f"{int(ANCHOR_P[b])}GPa": {"V": round(V, 4), "K": round(K, 2)}
                          for b, (V, K) in comps.items()},
           "fig1_digitized": {b: {"Fe": FIG1[b]["Fe"],
                                  "O": {str(k): v for k, v in FIG1[b]["O"].items()},
                                  "C": {str(k): v for k, v in FIG1[b]["C"].items()}} for b in FIG1},
           "source": "Badro, Cote & Brodholt 2014, PNAS 111, 7542, Fig. 1 (digitized)"}
    with open(os.path.join(HERE, "liquid_carbon_badro.json"), "w") as f:
        json.dump(out, f, indent=1)


if __name__ == "__main__":
    main()
