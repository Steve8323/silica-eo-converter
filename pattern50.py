#!/usr/bin/env python
"""The N=50 scaling pattern: per-tap coupler prescriptions.

Targets (exact ladder algebra, slots flat at 1/N^2):
    s_k = m_k = 1/(N+1-k),   k = 1..50   (k=50 -> full transfer)

Design rule (from the FDTD-calibrated anchors):
  * every coupler has a coupling FLOOR set by its adjacent arcs; floors and
    rates scale exponentially with gap (supermode + measured anchors):
        rate(g)      = (pi/2) / L_pi(g)            [supermode, exact]
        floor_tap(g) = 0.434 * rate(g)/rate(0.12)  [anchor: FDTD @0.12]
        floor_mrg(g) = 0.790 * rate(g)/rate(0.12)  [anchors: 0.79@0.12, 0.55@0.30 -> same shape]
  * per tap: choose the SMALLEST gap tier whose floor <= 0.85*kappa_target
    (15% trim headroom), then set the straight length
        L = (kappa_target - floor(g)) / rate(g),  clamped to [0.15, L_pi(g)]

Early taps (ratios down to 1/50) land in wide-gap tiers; late taps walk back
to gap 0.12 with growing L up to full transfer.  Output: fab table (JSON) +
pattern figure.  Final per-coupler values are then trimmed in situ exactly
as demonstrated on the N=6 device (2 passes -> +-1.5 dB).
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from multilayer import coupling_length
from dispersion import n_silica

N = 50
NG = 1.4948
C = 0.299792458
TIERS = [0.12, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]
rate = {g: np.pi / 2 / coupling_length(1.55, 0.991, g, n_silica(1.55))
        for g in TIERS}
r0 = rate[0.12]
FLOOR_T = {g: 0.434 * rate[g] / r0 for g in TIERS}
FLOOR_M = {g: 0.790 * rate[g] / r0 for g in TIERS}

def prescribe(kappa, floors):
    for g in TIERS:
        if floors[g] <= 0.85 * kappa:
            L = (kappa - floors[g]) / rate[g]
            Lpi = np.pi / 2 / rate[g]
            return g, float(np.clip(L, 0.15, Lpi))
    g = TIERS[-1]
    return g, 0.15

pitch = 14.0 + 0.3 * np.arange(N)
xs = 10.0 + np.concatenate([[0], np.cumsum(pitch[:-1])])
rows = []
for k in range(1, N + 1):
    ratio = 1.0 / (N + 1 - k)
    kap = float(np.arcsin(np.sqrt(min(ratio, 1.0))))
    gt, Lt = prescribe(kap, FLOOR_T)
    gm, Lm = prescribe(kap, FLOOR_M)
    rows.append(dict(k=k, ratio=ratio, kappa=kap,
                     tap_gap=gt, tap_L=round(Lt, 2),
                     mrg_gap=gm, mrg_L=round(Lm, 2),
                     x_um=round(float(xs[k - 1]), 1),
                     slot_fs=round((2 * (xs[k - 1] + 4) + 40) * NG / C, 0)))

with open("out/pattern50.json", "w") as f:
    json.dump(dict(note="s_k=m_k=1/(51-k); floors/rates from FDTD anchors; "
                        "in-situ trim (2 passes) finishes each coupler",
                   rate=rate, floor_tap=FLOOR_T, floor_mrg=FLOOR_M,
                   rows=rows), f, indent=1)

ks = [r["k"] for r in rows]
fig, axs = plt.subplots(3, 1, figsize=(9, 7.5), dpi=140, sharex=True)
axs[0].semilogy(ks, [r["ratio"] for r in rows], "o-", ms=3, color="#1f77b4")
axs[0].set_ylabel("target ratio s_k = m_k")
axs[0].set_title("N=50 ladder pattern: ratio -> gap tier -> length", fontsize=10)
axs[1].step(ks, [r["tap_gap"] for r in rows], where="mid", label="tap gap", color="#1f77b4")
axs[1].step(ks, [r["mrg_gap"] for r in rows], where="mid", label="merge gap", color="#d62728")
axs[1].set_ylabel("gap tier (um)"); axs[1].legend(fontsize=8)
axs[2].plot(ks, [r["tap_L"] for r in rows], "o-", ms=3, label="tap L", color="#1f77b4")
axs[2].plot(ks, [r["mrg_L"] for r in rows], "s-", ms=3, label="merge L", color="#d62728")
axs[2].set_ylabel("straight length (um)"); axs[2].set_xlabel("tap k")
axs[2].legend(fontsize=8)
fig.tight_layout(); fig.savefig("fig/pattern50.png"); plt.close(fig)

print(f"{'k':>3} {'ratio':>7} {'tapG':>5} {'tapL':>6} {'mrgG':>5} {'mrgL':>6} {'x':>7} {'slot_fs':>8}")
for r in rows[::7] + [rows[-1]]:
    print(f"{r['k']:3d} {r['ratio']:7.4f} {r['tap_gap']:5.2f} {r['tap_L']:6.2f} "
          f"{r['mrg_gap']:5.2f} {r['mrg_L']:6.2f} {r['x_um']:7.1f} {r['slot_fs']:8.0f}")
