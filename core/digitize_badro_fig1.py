"""Digitize the C and O markers of Badro et al. (2014) Fig. 1 (provenance for liquid_carbon_badro.py).

Usage (needs PyMuPDF; it is in environment.yml):
    python digitize_badro_fig1.py <Badro_2014 main-text PDF>
Prints colour-blob centroids converted to data values. Axis calibration comes from the tick marks
of the 895 x 846 px Fig. 1 raster on page 2. Overlapping or occluded markers (O under C at ICB
density; C at 8.3 mol% under the CMB velocity target line) were resolved by inspecting enlarged
crops, as recorded in liquid_carbon_badro.py.
"""
import io
import sys

import fitz
import json
import numpy as np
from PIL import Image
doc = fitz.open(sys.argv[1])
raw = doc.extract_image(doc[1].get_images(full=True)[1][0])["image"]
a = np.asarray(Image.open(io.BytesIO(raw)).convert("RGB")).astype(int)
R, G, B = a[..., 0], a[..., 1], a[..., 2]
masks = {
    "C": (abs(R - G) < 18) & (abs(G - B) < 18) & (R > 90) & (R < 175),   # grey hexagons
    "O": (R > 170) & (G < 90) & (B < 90),                                  # red circles
    "any": (a.sum(2) < 600),                                               # any marker colour at x = 0
}
# panel: (x0 pixel at 0 mol%, px per mol%, y_top pixel, value_top, y_bottom pixel, value_bottom)
panels = {
    "CMB_rho": (85, 68 / 6, 7, 10800, 357, 9200),
    "CMB_v":   (546, 68 / 6, 7, 8500, 357, 7700),
    "ICB_rho": (85, 68 / 6, 446, 13400, 796, 11400),
    "ICB_v":   (546, 68 / 6, 446, 10700, 796, 9900),
}
def yval(y, p):
    _, _, yt, vt, yb, vb = p
    return vt + (y - yt) * (vb - vt) / (yb - yt)
out = {}
for name, p in panels.items():
    x0, ppm = p[0], p[1]
    for el in ("C", "O"):
        for conc in (8.3, 16.7):
            xc = int(round(x0 + ppm * conc))
            m = masks[el][:, xc - 6:xc + 7]
            ys = [y for y in range(p[2], p[4]) if m[y].sum() >= 6]   # wide rows = marker body
            if not ys:
                out[f"{name}_{el}_{conc}"] = None; continue
            # largest contiguous group
            groups = [[ys[0]]]
            for y in ys[1:]:
                (groups[-1].append(y) if y - groups[-1][-1] <= 1 else groups.append([y]))
            g = max(groups, key=len)
            out[f"{name}_{el}_{conc}"] = dict(y_px=float(np.mean(g)), value=round(yval(np.mean(g), p), 1), rows=len(g))
    m = masks["any"][:, x0 - 6:x0 + 7]
    ys = [y for y in range(p[2] + 3, p[4] - 3) if m[y].sum() >= 7]
    groups = [[ys[0]]]
    for y in ys[1:]:
        (groups[-1].append(y) if y - groups[-1][-1] <= 1 else groups.append([y]))
    g = max(groups, key=len)
    out[f"{name}_Fe_0"] = dict(y_px=float(np.mean(g)), value=round(yval(np.mean(g), p), 1), rows=len(g))
print(json.dumps(out, indent=1))
