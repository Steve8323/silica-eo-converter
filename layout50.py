#!/usr/bin/env python
"""Full to-scale circuit layout, N=50 silica-only EO->photonic converter,
with electrodes and optical ports labeled.

Signal path:  CW laser (left facet) -> master MRR (fast EO electrode pair,
GHz drive + DC bias) -> zero-GVD bus w=0.991 -> 50-tap counter-directional
ladder (weight electrode pair on every arm) -> serialized output (left facet
of bus B).

Electrodes are coplanar pairs (signal + ground) flanking the silica guide;
they act via the DC Kerr effect (eo_mapping.py).  Dumps all coordinates +
coupler prescriptions to out/layout50.json.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle, FancyArrow

N = 50
W = 0.991
RHO = 6.0
MRR_R = 16.0
X_MRR = 60.0
X_LADDER = 260.0

with open("out/pattern50.json") as f:
    pat = json.load(f)["rows"]
xs = np.array([r["x_um"] for r in pat])
L1 = np.array([r["tap_L"] for r in pat])

yA = 0.0
y_arm = yA - 1.2
y_mrg = y_arm - 2 * RHO
yB = y_mrg - 1.4
x_end = X_LADDER + xs[-1] + L1[-1] + 2 * RHO + 10

GOLD = "#b9bec7"  # silver
fig, ax = plt.subplots(figsize=(14, 5.4), dpi=150)

# ---- optical ports ----
ax.add_patch(FancyArrow(-46, yA, 30, 0, width=1.6, head_width=5, head_length=7,
                        color="#c0392b", alpha=0.9))
ax.text(-30, yA + 9, "CW LASER IN\n(1550 nm facet)", ha="center", fontsize=7.5,
        color="#c0392b", weight="bold")
ax.add_patch(FancyArrow(X_LADDER - 60, yB, -30, 0, width=1.6, head_width=5,
                        head_length=7, color="#2ca02c", alpha=0.9))
ax.text(85, yB - 20, "SERIALIZED OUTPUT\n(50 slots, 0.3-10.7 ps, then blank)",
        ha="center", fontsize=7.5, color="#2ca02c", weight="bold")

# ---- buses ----
ax.plot([0, x_end], [yA, yA], color="#1f77b4", lw=1.6)
ax.plot([X_LADDER - 60, x_end - 8], [yB, yB], color="#2ca02c", lw=1.6)
ax.text(X_LADDER + 20, yA + 3, "bus A: input pulse ->", fontsize=7, color="#1f77b4")
ax.text(X_LADDER + 150, yB - 7, "<- bus B: output", fontsize=7, color="#2ca02c")
ax.text((X_MRR + X_LADDER) / 2 + 30, yA + 4, "zero-GVD bus w=0.991 um",
        fontsize=7, ha="center")

# ---- master MRR + fast EO electrodes ----
cy = yA + W / 2 + 0.15 + MRR_R
ax.add_patch(Circle((X_MRR, cy), MRR_R, fill=False, lw=1.4, color="#1f77b4"))
ax.add_patch(Circle((X_MRR, cy), MRR_R + 3.5, fill=False, lw=4.0, color=GOLD, alpha=0.8))
ax.add_patch(Circle((X_MRR, cy), MRR_R - 3.5, fill=False, lw=4.0, color=GOLD, alpha=0.8))
ax.add_patch(Rectangle((X_MRR + MRR_R + 2, cy + 8), 14, 10, color=GOLD, alpha=0.9))
ax.add_patch(Rectangle((X_MRR - MRR_R - 16, cy - 18), 14, 10, color=GOLD, alpha=0.9))
ax.text(X_MRR + MRR_R + 22, cy + 13, "RF pad (GHz drive)", ha="left", fontsize=6.5)
ax.text(X_MRR - MRR_R - 20, cy - 32, "DC bias pad\n(E_b linearizes Kerr)", ha="right", fontsize=6.5)
ax.text(X_MRR + 150, cy + MRR_R + 14,
        "MASTER MRR R=16 um (Q_bend=1.1e7, FSR 2 THz)\n"
        "FAST EO ELECTRODES: inner+outer Ag rings flanking the guide\n"
        "DC-Kerr tuning, ~0.1-1 GHz photon-lifetime bandwidth",
        ha="center", fontsize=7, weight="bold")

# ---- ladder cells + weight electrodes ----
for k in range(N):
    x = X_LADDER + xs[k]
    xu = x + L1[k]
    ax.plot([x, xu], [y_arm, y_arm], color="#444", lw=0.9)
    th = np.linspace(-np.pi / 2, np.pi / 2, 20)
    ax.plot(xu + RHO * np.cos(th), y_arm - RHO + RHO * np.sin(th), color="#444", lw=0.9)
    ax.plot([xu - max(L1[k], 2), xu], [y_mrg, y_mrg], color="#444", lw=0.9)
    ax.plot([xu - max(L1[k], 2) - 4, xu - max(L1[k], 2) - 4], [y_mrg - 4, yB - 8],
            color="#bbb", lw=0.6)
    # weight electrode PAIR on the arm's vertical run (right side of U)
    ax.add_patch(Rectangle((xu + RHO + 0.6, y_mrg + 2), 1.4, 2 * RHO - 4,
                           color=GOLD, alpha=0.85, lw=0))
    ax.add_patch(Rectangle((xu + RHO - 2.9, y_mrg + 2), 1.4, 2 * RHO - 4,
                           color=GOLD, alpha=0.45, lw=0))
kk = 24
xk = X_LADDER + xs[kk] + L1[kk] + RHO
ax.annotate("Ag WEIGHT ELECTRODE PAIR on every arm\n(signal + ground; slow DC-Kerr "
            "phase/amplitude trim of slot k)",
            xy=(xk, y_mrg + RHO), xytext=(xk - 120, y_mrg - 26),
            fontsize=7, weight="bold", color="#8a6d00",
            arrowprops=dict(arrowstyle="->", lw=0.8, color="#8a6d00"))
ax.text(X_LADDER + xs[-1] / 2, yA + 10,
        "50-TAP LADDER: tap ratios 1/50 -> 1 (pattern50.json), uneven pitch "
        "14 -> 28.7 um -> uneven slot times", ha="center", fontsize=7)

# ---- zoom inset of one cell ----
axi = fig.add_axes([0.72, 0.56, 0.155, 0.36])
k = 30
x, l1 = 0.0, 3.0
axi.plot([-4, 10], [0, 0], color="#1f77b4", lw=2)
axi.plot([-4, 10], [-14, -14], color="#2ca02c", lw=2)
axi.plot([x - 1, x + l1], [-1.2, -1.2], color="#444", lw=1.5)
th = np.linspace(-np.pi / 2, np.pi / 2, 30)
axi.plot(x + l1 + 5 * np.cos(th), -6.2 + 5 * np.sin(th), color="#444", lw=1.5)
axi.plot([x - 2, x + l1], [-11.2, -11.2], color="#444", lw=1.5)
axi.plot([x - 4, x - 4], [-13.5, -20], color="#bbb", lw=1)
axi.add_patch(Rectangle((x + l1 + 5.8, -10), 1.0, 7.5, color=GOLD, alpha=0.9))
axi.add_patch(Rectangle((x + l1 + 2.4, -10), 1.0, 7.5, color=GOLD, alpha=0.45))
axi.text(1.5, 1.2, "tap coupler", fontsize=6)
axi.text(-3.6, -10.6, "merge coupler", fontsize=6)
axi.text(x + l1 + 7.2, -6.5, "weight\nelectrodes", fontsize=6, color="#8a6d00")
axi.text(-3.8, -16.8, "dump", fontsize=6, color="#999")
axi.set_title("one ladder cell (zoom)", fontsize=7)
axi.set_xticks([]); axi.set_yticks([])
axi.set_xlim(-6, 14); axi.set_ylim(-21, 4); axi.set_aspect("equal")

ax.set_xlim(-60, x_end + 15)
ax.set_ylim(yB - 36, cy + MRR_R + 44)
ax.set_aspect("equal")
ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
ax.set_title("Silica-only EO->photonic converter, N=50 — layout with electrodes "
             "and optical ports labeled (to scale)", fontsize=9, loc="left")
fig.subplots_adjust(left=0.05, right=0.99, top=0.97, bottom=0.12)
fig.savefig("fig/layout50.png"); plt.close(fig)

table = [dict(r, x_um=float(X_MRR if False else X_LADDER + xs[i])) for i, r in enumerate(pat)]
with open("out/layout50.json", "w") as f:
    json.dump(dict(W=W, rho=RHO, mrr_R=MRR_R, mrr_x=X_MRR,
                   laser_in="left facet, bus A, 1550 nm CW",
                   output="left facet, bus B (serialized slots)",
                   electrodes=dict(
                       mrr_fast="coplanar ring pair flanking master MRR, RF pad + DC bias pad",
                       weights="one signal+ground pair per ladder arm, slow DC-Kerr trim"),
                   cells=table), f, indent=1)
print("layout50.png + layout50.json updated with labeled electrodes/ports")
