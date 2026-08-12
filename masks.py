#!/usr/bin/env python
"""Photomask artwork for the N=50 silica-only converter.

Three mask layers, all coordinates in um, generated from the calibrated
pattern50.json prescriptions:

  WG  (mask SEO-WG-01, dark-field, etch): waveguides — buses, master MRR,
      50 ladder cells (per-cell tap gap/length, U-bend, merge gap/length,
      dump), facet tapers, alignment marks.
  MET (mask SEO-MET-02, clear-field, liftoff): master-MRR ring electrode
      pair + RF/DC pads, per-arm weight electrode pairs + routing rails.
  REL (mask SEO-REL-03, dark-field, release): XeF2 access windows to
      undercut the Si and suspend the silica (air cladding on both sides —
      the 'silica only' requirement).

Outputs: out/masks.json (feature lists per layer for the viewer),
fig/mask_*.png previews.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

W = 0.991
RHO_MIN = 6.0
MRR_R = 16.0
MRR_GAP = 0.15
X_MRR = 60.0
X_LADDER = 260.0
X_FACET = -120.0
E_W = 3.0          # electrode width
E_GAP = 1.0        # electrode edge to waveguide edge

with open("out/pattern50.json") as f:
    rows = json.load(f)["rows"]

WGf, METf, RELf = [], [], []   # feature lists

def line(layer, pts, w=W):
    layer.append(dict(t="line", p=[[round(a, 3), round(b, 3)] for a, b in pts], w=w))

def arc(layer, c, r, a0, a1, w=W):
    layer.append(dict(t="arc", c=[round(c[0], 3), round(c[1], 3)], r=round(r, 3),
                      a0=a0, a1=a1, w=w))

def circ(layer, c, r, w):
    arc(layer, c, r, 0, 360, w)

def rect(layer, x, y, dx, dy):
    layer.append(dict(t="rect", x=round(x, 3), y=round(y, 3),
                      dx=round(dx, 3), dy=round(dy, 3)))

# ---------------- WG layer ----------------
yA = 0.0
yB = -(W + 0.70) - 2 * RHO_MIN - (W + 0.70) - 1.0   # deepest cell fits
x_end = X_LADDER + rows[-1]["x_um"] - 10.0 + rows[-1]["tap_L"] + 2 * RHO_MIN + 20
# facet tapers (inverse taper 0.3 -> W over 60 um)
line(WGf, [(X_FACET, yA), (X_FACET + 60, yA)], w=0.3)
line(WGf, [(X_FACET + 60, yA), (X_FACET + 120, yA)], w=0.65)
line(WGf, [(X_FACET + 120, yA), (x_end, yA)])                       # bus A
line(WGf, [(X_FACET + 120, yB), (x_end - 10, yB)])                  # bus B
line(WGf, [(X_FACET, yB), (X_FACET + 60, yB)], w=0.3)
line(WGf, [(X_FACET + 60, yB), (X_FACET + 120, yB)], w=0.65)
# master MRR (ring above bus A)
mrr_c = (X_MRR, yA + W / 2 + MRR_GAP + MRR_R + W / 2)
circ(WGf, mrr_c, MRR_R, W)

cells = []
for r in rows:
    x = X_LADDER + r["x_um"] - 10.0     # pattern50 x is ladder-relative (starts at 10)
    l1, l2 = r["tap_L"], r["mrg_L"]
    g1, g2 = r["tap_gap"], r["mrg_gap"]
    y_arm = yA - (W + g1)
    y_mrg = yB + (W + g2)
    rho = (y_arm - y_mrg) / 2
    xu = x + l1
    cy = y_arm - rho
    line(WGf, [(x - 1.5, y_arm), (xu, y_arm)])                      # tap section
    arc(WGf, (xu, cy), rho, -90, 90)                                # U-bend
    line(WGf, [(xu - max(l2, 1.5), y_mrg), (xu, y_mrg)])            # merge section
    xd = xu - max(l2, 1.5) - 3.0
    line(WGf, [(xd, y_mrg), (xu - max(l2, 1.5), y_mrg)])            # lead
    arc(WGf, (xd, y_mrg - 4.0), 4.0, 90, 180)                       # down bend
    line(WGf, [(xd - 4.0, y_mrg - 4.0), (xd - 4.0, yB - 10)])       # dump stub
    cells.append(dict(x=x, xu=xu, y_arm=y_arm, y_mrg=y_mrg, rho=rho, cy=cy))

# ---------------- MET layer ----------------
circ(METf, mrr_c, MRR_R + W / 2 + E_GAP + E_W / 2, E_W)             # outer ring el.
circ(METf, mrr_c, MRR_R - W / 2 - E_GAP - E_W / 2, E_W)             # inner ring el.
rect(METf, mrr_c[0] + MRR_R + 6, mrr_c[1] + 6, 20, 16)              # RF pad
line(METf, [(mrr_c[0] + MRR_R + E_GAP + E_W / 2, mrr_c[1] + 4),
            (mrr_c[0] + MRR_R + 10, mrr_c[1] + 10)], w=2)
rect(METf, mrr_c[0] - MRR_R - 26, mrr_c[1] + 6, 20, 16)             # DC pad
line(METf, [(mrr_c[0] - MRR_R + E_GAP - E_W / 2 + 1, mrr_c[1] + 2),
            (mrr_c[0] - MRR_R - 12, mrr_c[1] + 10)], w=2)
# weight electrodes: pair flanking each arm's vertical run + rails
rail_y = 16.0
line(METf, [(X_LADDER - 20, rail_y), (x_end, rail_y)], w=6)         # signal rail
line(METf, [(X_LADDER - 20, rail_y + 12), (x_end, rail_y + 12)], w=6)  # ground rail
rect(METf, X_LADDER - 60, rail_y - 4, 30, 30)                       # rail pad
for i, c in enumerate(cells):
    xs_el = c["xu"] + c["rho"] + W / 2 + E_GAP
    yl, yh = c["y_mrg"] + 2, c["y_arm"] - 2
    rect(METf, xs_el, yl, E_W, yh - yl)                             # signal finger
    rect(METf, xs_el - W - 2 * E_GAP - E_W - 2 * c["rho"], yl, E_W, yh - yl)  # gnd finger
    line(METf, [(xs_el + E_W / 2, yh), (xs_el + E_W / 2, rail_y)], w=1.2)
    line(METf, [(xs_el - W - 2 * E_GAP - E_W / 2 - 2 * c["rho"], yl),
                (xs_el - W - 2 * E_GAP - E_W / 2 - 2 * c["rho"], yB - 14),
                (x_end + 6, yB - 14)], w=1.2)

# ---------------- REL layer (release windows) ----------------
for x in np.arange(X_FACET + 130, x_end - 5, 24):
    rect(RELf, x, yA + 3.0, 10, 6)                                  # above bus A
    rect(RELf, x, yB - 9.0, 10, 6)                                  # below bus B
for c in cells[::2]:
    rect(RELf, c["x"] + 2, c["cy"] - 3, min(2 * c["rho"] - 10, 8), 6)  # inside U

# alignment marks (all layers)
for L in (WGf, METf, RELf):
    for mx, my in ((X_FACET + 10, 60), (x_end - 10, 60), (X_FACET + 10, -80), (x_end - 10, -80)):
        line(L, [(mx - 8, my), (mx + 8, my)], w=1.5)
        line(L, [(mx, my - 8), (mx, my + 8)], w=1.5)

masks = dict(
    units="um",
    chip=dict(x0=X_FACET - 20, x1=x_end + 40, y0=-100, y1=80),
    layers=[
        dict(name="WG  (SEO-WG-01, etch)", color="#1f77b4", features=WGf),
        dict(name="MET (SEO-MET-02, liftoff)", color="#d4a017", features=METf),
        dict(name="REL (SEO-REL-03, release)", color="#e377c2", features=RELf),
    ])
with open("out/masks.json", "w") as f:
    json.dump(masks, f)

# ---- PNG previews (one per mask, dark-field style) ----
def draw(ax, feats, color):
    for ft in feats:
        if ft["t"] == "line":
            p = np.array(ft["p"])
            ax.plot(p[:, 0], p[:, 1], color=color, lw=0.8, solid_capstyle="butt")
        elif ft["t"] == "arc":
            th = np.radians(np.linspace(ft["a0"], ft["a1"], 60))
            ax.plot(ft["c"][0] + ft["r"] * np.cos(th),
                    ft["c"][1] + ft["r"] * np.sin(th), color=color, lw=0.8)
        else:
            ax.add_patch(plt.Rectangle((ft["x"], ft["y"]), ft["dx"], ft["dy"],
                                       color=color, alpha=0.8, lw=0))

for lay in masks["layers"]:
    fig, ax = plt.subplots(figsize=(14, 3.2), dpi=150)
    ax.set_facecolor("black")
    draw(ax, lay["features"], lay["color"])
    ax.set_xlim(masks["chip"]["x0"], masks["chip"]["x1"])
    ax.set_ylim(masks["chip"]["y0"], masks["chip"]["y1"])
    ax.set_aspect("equal")
    ax.set_title(f"photomask {lay['name']}", fontsize=9)
    ax.set_xlabel("x (um)")
    tag = lay["name"].split()[0].lower()
    fig.tight_layout(); fig.savefig(f"fig/mask_{tag}.png"); plt.close(fig)
print("masks.json +", ", ".join(f"mask_{l['name'].split()[0].lower()}.png"
                                for l in masks["layers"]))
