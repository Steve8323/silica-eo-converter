#!/usr/bin/env python
"""Figures for the all-pass ring-chain serializer runs."""
import glob
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def J(p):
    with open(p) as f:
        return json.load(f)

files = sorted(glob.glob("out/serializer_N*_g0.1.json")) or \
        sorted(glob.glob("out/serializer_N*.json"))
runs = sorted([J(f) for f in files], key=lambda d: d["N"])
if not runs:
    raise SystemExit("no serializer runs")

fig, axs = plt.subplots(len(runs), 1, figsize=(9, 1.9 * len(runs) + 1),
                        dpi=140, sharex=True)
axs = np.atleast_1d(axs)
for ax, d in zip(axs, runs):
    t = np.array(d["t_fs"]); e = np.array(d["env"])
    e = e / e.max()
    ax.semilogy(t, np.maximum(e, 1e-5), color="#1f77b4", lw=0.9)
    ax.fill_between(t, np.maximum(e, 1e-5), 1e-5, color="#1f77b4", alpha=0.5, lw=0)
    ax.set_ylabel(f"N={d['N']}\n(norm., log)", fontsize=8)
    ax.set_ylim(1e-4, 1.5)
    ax.text(0.99, 0.85, f"E_out/E_in = {d['ratio']:.3f}",
            transform=ax.transAxes, ha="right", fontsize=8,
            color="#d62728")
axs[-1].set_xlabel("time at output port (fs)")
axs[0].set_title("All-pass ring chain: uneven pulse cluster then blank "
                 "(MEEP time-domain output)", fontsize=10)
fig.tight_layout(); fig.savefig("fig/serializer_trains.png"); plt.close(fig)

fig, ax = plt.subplots(figsize=(5.5, 3.4), dpi=140)
Ns = [d["N"] for d in runs]; rat = [d["ratio"] for d in runs]
ax.plot(Ns, rat, "o-", color="#1f77b4")
ax.axhline(1.0, color="k", lw=0.6, ls="--")
ax.set_ylim(0, 1.1)
ax.set_xlabel("number of rings N"); ax.set_ylabel("E_out / E_in")
ax.set_title("Zero-loss scaling check (all-pass chain)", fontsize=10)
ax.set_xticks(Ns)
fig.tight_layout(); fig.savefig("fig/serializer_scaling.png"); plt.close(fig)
print("saved serializer figs;", [(d["N"], round(d["ratio"], 4)) for d in runs])
