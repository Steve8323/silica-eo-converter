#!/usr/bin/env python
"""N=50 tap-ladder: exact FIR composition + predicted output waveform.

Slot transmission (exact ladder algebra, telescoping):
  T_k = [s_k prod_{j<k}(1-s_j)] * [m_k prod_{j<k}(1-m_j)] = 1/N^2  (uniform)
with s_k = m_k = 1/(N-k), k=0..N-1.  On top of the fundamental 1/N^2 we
apply MEASURED per-element excess: through-coupler excess and the
perpendicular dump-stub crossings each slot passes.  FDTD anchor: the N=6
run (serializer_ladder.py) validates slot positions/uniformity; this script
composes the same algebra to N=50.
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N = 50
NG = 1.4948
C_UMFS = 0.299792458
T_CPL = 0.998        # through-coupler excess transmission (per pass)
T_XNG = 0.995        # perpendicular crossing transmission (2D FDTD-level)
PULSE_FS = 35.0

pitch = 14.0 + 0.3 * np.arange(N)          # uneven, mildly chirped (um)
xs = 10.0 + np.concatenate([[0], np.cumsum(pitch[:-1])])
slot_t = 2 * (xs + 4) * NG / C_UMFS        # fs
s = np.array([1.0 / (N - k) for k in range(N)])
T_fund = np.zeros(N)
for k in range(N):
    split = s[k] * np.prod(1 - s[:k])
    merge = s[k] * np.prod(1 - s[:k])
    T_fund[k] = split * merge
excess = np.array([T_CPL ** (2 * k) * T_XNG ** k for k in range(N)])
T = T_fund * excess

t = np.arange(0, slot_t[-1] + 3000, 4.0)
env = np.zeros_like(t)
for k in range(N):
    env += T[k] * np.exp(-((t - slot_t[k]) / (PULSE_FS / 1.665)) ** 2)

fig, axs = plt.subplots(2, 1, figsize=(10, 5.6), dpi=140,
                        gridspec_kw=dict(height_ratios=[2, 1]))
axs[0].semilogy(t / 1000, np.maximum(env / env.max(), 1e-4), color="#1f77b4", lw=0.8)
axs[0].fill_between(t / 1000, np.maximum(env / env.max(), 1e-4), 1e-4,
                    color="#1f77b4", alpha=0.5)
axs[0].set_xlabel("time (ps)"); axs[0].set_ylabel("output power (norm., log)")
axs[0].set_title(f"N=50 ladder output: 50 sequential slots then EXACTLY blank "
                 f"(FIR — no recirculation)", fontsize=10)
axs[0].axvspan(slot_t[-1] / 1000 + 0.15, t[-1] / 1000, color="green", alpha=0.06)
axs[0].text(0.99, 0.9, "blank tail (structural)", transform=axs[0].transAxes,
            ha="right", fontsize=8, color="green")
axs[1].plot(np.arange(1, N + 1), T / T.mean(), "o-", ms=3, color="#d62728")
axs[1].set_xlabel("slot k"); axs[1].set_ylabel("T_k / mean")
axs[1].set_title("slot uniformity: fundamental exactly 1/N^2; residual slope = "
                 "excess (couplers+crossings) accumulating toward late slots",
                 fontsize=9)
fig.tight_layout(); fig.savefig("fig/ladder50_waveform.png"); plt.close(fig)

res = dict(N=N, T_fund_slot=float(T_fund[0]), total_fund=float(T_fund.sum()),
           total_with_excess=float(T.sum()),
           slot_first_fs=float(slot_t[0]), slot_last_fs=float(slot_t[-1]),
           uniformity_max_over_min=float(T.max() / T.min()))
with open("out/ladder50_model.json", "w") as f:
    json.dump(dict(res=res, slot_t=list(slot_t), T=list(T)), f)
print(f"N=50: per-slot fundamental T = {T_fund[0]:.6f} (=1/2500), "
      f"total = {T_fund.sum():.4f} (=1/50)")
print(f"with excess: total {T.sum():.4f}; slot spread {res['uniformity_max_over_min']:.2f}x")
print(f"train: {slot_t[0]/1000:.2f} -> {slot_t[-1]/1000:.2f} ps, then blank")
