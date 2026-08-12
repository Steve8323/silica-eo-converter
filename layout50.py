#!/usr/bin/env python
"""Full to-scale circuit layout, N=50 silica-only EO->photonic converter.

  [CW in] -> master MRR (R=16, fast EO electrode) -> zero-GVD bus (w=0.991)
  -> 50-tap counter-directional ladder (split ladder on bus A, mirrored
  merge ladder on bus B, per-arm weight electrode) -> serialized output.

Draws every waveguide to scale and dumps all coordinates + coupler lengths
to out/layout50.json (the fab table).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Rectangle

N = 50
W = 0.991
GAP = 0.12          # tap side; merge side uses 0.40 (arc-coupling floor)
LPI = 11.31
RHO = 6.0

def clen(r):
    return LPI * 2 / np.pi * np.arcsin(np.sqrt(min(r, 1.0)))

pitch = 14.0 + 0.3 * np.arange(N)
xs = np.concatenate([[0], np.cumsum(pitch[:-1])])
s = np.array([1.0 / (N - k) for k in range(N)])
L1 = np.array([clen(r) for r in s])

MRR_R = 16.0
X_MRR = 60.0
X_LADDER = 260.0           # ladder start (after mode settling)
yA, y_arm = 0.0, -(W + GAP)
y_mrg = y_arm - 2 * RHO
yB = y_mrg - (W + GAP)

fig, ax = plt.subplots(figsize=(13, 3.6), dpi=150)
# master MRR + its electrode
ax.add_patch(Circle((X_MRR, yA + W / 2 + 0.15 + MRR_R), MRR_R, fill=False,
                    lw=1.2, color="#1f77b4"))
ax.add_patch(Circle((X_MRR, yA + W / 2 + 0.15 + MRR_R), MRR_R + 4, fill=False,
                    lw=3.0, color="#f0b000", alpha=0.5))
ax.text(X_MRR, yA + MRR_R + 24, "master MRR R=16 um\nQ_bend=1.1e7, FSR 2 THz\n"
        "fast EO electrode (GHz)", ha="center", fontsize=6.5)
# buses
x_end = X_LADDER + xs[-1] + L1[-1] + 2 * RHO + 10
ax.plot([0, x_end], [yA, yA], color="#1f77b4", lw=1.4)
ax.plot([0, x_end - 8], [yB, yB], color="#2ca02c", lw=1.4)
ax.annotate("bus A (input pulse ->)", (X_LADDER - 90, yA + 2), fontsize=7)
ax.annotate("<- bus B (serialized output)", (X_LADDER - 110, yB - 4.5),
            fontsize=7, color="#2ca02c")
# ladder cells
for k in range(N):
    x = X_LADDER + xs[k]
    xu = x + L1[k]
    ax.plot([x, xu], [y_arm, y_arm], color="#444", lw=1.0)
    th = np.linspace(-np.pi / 2, np.pi / 2, 24)
    ax.plot(xu + RHO * np.cos(th), y_arm - RHO + RHO * np.sin(th),
            color="#444", lw=1.0)
    ax.plot([xu - L1[k], xu], [y_mrg, y_mrg], color="#444", lw=1.0)
    xd = xu - L1[k] - 2.0
    ax.plot([xd, xu - L1[k]], [y_mrg, y_mrg], color="#444", lw=0.8)
    ax.plot([xd - 2, xd - 2], [y_mrg - 2, yB - 6], color="#999", lw=0.7)
    # weight electrode on the arm (phase/amplitude trim)
    ax.add_patch(Rectangle((x + 1, y_arm + 0.8), max(L1[k] - 2, 3), 1.2,
                           color="#f0b000", alpha=0.6, lw=0))
ax.text(X_LADDER + xs[-1] / 2, y_arm + 8,
        f"50-tap ladder: split s_k=1/(50-k), mirrored merge; "
        f"uneven pitch 14->28.7 um; weight electrode per arm (gold)",
        ha="center", fontsize=7)
ax.annotate("", xy=(X_LADDER, yA + 6), xytext=(X_MRR + MRR_R + 8, yA + 6),
            arrowprops=dict(arrowstyle="->", lw=0.8))
ax.text((X_LADDER + X_MRR) / 2 + 20, yA + 8, "zero-GVD bus w=0.991 um",
        fontsize=6.5, ha="center")
ax.set_xlim(-15, x_end + 15)
ax.set_ylim(yB - 16, yA + MRR_R * 2 + 22)
ax.set_aspect("equal")
ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
ax.set_title(f"Silica-only EO->photonic converter, N=50 — full layout to scale "
             f"({x_end:.0f} x {MRR_R*2+38+abs(yB):.0f} um core; "
             f"5 mm chip budget: {x_end/5000*100:.0f}% used)", fontsize=9)
fig.tight_layout(); fig.savefig("fig/layout50.png"); plt.close(fig)

table = [dict(k=k + 1, x_um=float(X_LADDER + xs[k]), s_ratio=float(s[k]),
              L_coupler_um=float(L1[k]),
              slot_fs=float(2 * (xs[k] + 4) * 1.4948 / 0.299792458))
         for k in range(N)]
with open("out/layout50.json", "w") as f:
    json.dump(dict(W=W, gap=GAP, Lpi=LPI, rho=RHO, mrr_R=MRR_R,
                   yA=yA, yB=yB, cells=table), f, indent=1)
print(f"layout: {x_end:.0f} um long; couplers L from {L1[0]:.2f} to {L1[-1]:.2f} um;"
      f" table -> out/layout50.json")
