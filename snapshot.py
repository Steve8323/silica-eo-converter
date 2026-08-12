#!/usr/bin/env python
"""Field snapshots: what the FDTD actually looks like.

1) MZI at CW 1550 with DN=0 (all power to cross port) and DN=0.0332 (bar port)
2) Master MRR whispering-gallery mode (R=8) ringing after pulsed excitation

usage: python snapshot.py mzi DN | ring
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import meep as mp

MODE = sys.argv[1]
W = 0.991
NSIL = 1.444024
RES = 20

if MODE == "mzi":
    DN = float(sys.argv[2])
    GAP_C, LC, LB, LP = 0.30, 5.35, 14.0, 24.0
    D_C, D_A, PADX, LIN = W + GAP_C, 6.0, 4.0, 3.0
    seg = [("in", LIN), ("b1", LB), ("c1", LC), ("b2", LB), ("ph", LP),
           ("b3", LB), ("c2", LC), ("b4", LB), ("out", LIN)]
    edges = np.cumsum([0] + [s[1] for s in seg])
    LTOT = edges[-1]
    segx = {n: (edges[i], edges[i + 1]) for i, (n, _) in enumerate(seg)}

    def y_top(x):
        for name, (a, b) in segx.items():
            if a <= x <= b:
                t = (x - a) / (b - a)
                c = (1 - np.cos(np.pi * t)) / 2
                if name in ("in", "ph", "out"):
                    return D_A / 2
                if name in ("b1", "b3"):
                    return D_A / 2 - (D_A - D_C) / 2 * c
                if name in ("c1", "c2"):
                    return D_C / 2
                return D_C / 2 + (D_A - D_C) / 2 * c
        return D_A / 2

    def material(p):
        x = min(max(p.x + LTOT / 2, 0.0), LTOT)
        yc = y_top(x)
        if abs(p.y - yc) <= W / 2:
            dn = DN if segx["ph"][0] <= x <= segx["ph"][1] else 0.0
            return mp.Medium(index=NSIL + dn)
        if abs(p.y + yc) <= W / 2:
            return mp.Medium(index=NSIL)
        return mp.Medium(index=1.0)

    sx, sy = LTOT + 2 * PADX, D_A + W + 10
    sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), resolution=RES,
                        boundary_layers=[mp.PML(1.5)],
                        material_function=material,
                        extra_materials=[mp.Medium(index=NSIL + DN + 1e-4)],
                        default_material=mp.Medium(index=1.0))
    sim.sources = [mp.EigenModeSource(
        mp.ContinuousSource(1 / 1.55, width=20),
        center=mp.Vector3(-sx / 2 + 2.2, D_A / 2), size=mp.Vector3(0, 3.6),
        eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]
    sim.run(until=400)
    eps = sim.get_array(component=mp.Dielectric)
    ez = sim.get_array(component=mp.Ez)
    fig, ax = plt.subplots(figsize=(12, 2.8), dpi=150)
    ax.imshow(eps.T, origin="lower", cmap="gray_r", alpha=0.35,
              extent=[-sx / 2, sx / 2, -sy / 2, sy / 2], aspect="auto")
    m = np.abs(ez).max()
    ax.imshow(ez.T, origin="lower", cmap="RdBu", vmin=-m, vmax=m, alpha=0.8,
              extent=[-sx / 2, sx / 2, -sy / 2, sy / 2], aspect="auto")
    ax.set_title(f"MZI steady state, electrode dn={DN:g} "
                 f"({'cross port (s=-1)' if DN == 0 else 'bar port (s=+1)'})",
                 fontsize=10)
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    fig.tight_layout()
    fig.savefig(f"fig/field_mzi_dn{DN:g}.png")
    print("saved", f"fig/field_mzi_dn{DN:g}.png")

else:  # ring
    R = 8.0
    pad = 3.0
    s = 2 * (R + W + pad)
    geom = [mp.Cylinder(radius=R + W / 2, material=mp.Medium(index=NSIL)),
            mp.Cylinder(radius=R - W / 2, material=mp.Medium(index=1.0))]
    sim = mp.Simulation(cell_size=mp.Vector3(s, s), geometry=geom,
                        resolution=RES, boundary_layers=[mp.PML(1.5)])
    sim.sources = [mp.Source(mp.GaussianSource(1 / 1.5486, fwidth=0.002),
                             component=mp.Ez, center=mp.Vector3(R, 0.11))]
    sim.run(until=800)
    eps = sim.get_array(component=mp.Dielectric)
    ez = sim.get_array(component=mp.Ez)
    fig, ax = plt.subplots(figsize=(5.5, 5.2), dpi=150)
    ax.imshow(eps.T, origin="lower", cmap="gray_r", alpha=0.35,
              extent=[-s / 2, s / 2, -s / 2, s / 2])
    m = np.abs(ez).max()
    ax.imshow(ez.T, origin="lower", cmap="RdBu", vmin=-m, vmax=m, alpha=0.8,
              extent=[-s / 2, s / 2, -s / 2, s / 2])
    ax.set_title("Master MRR: whispering-gallery mode, R=8 um (1548.6 nm)",
                 fontsize=10)
    ax.set_xlabel("x (um)"); ax.set_ylabel("y (um)")
    fig.tight_layout()
    fig.savefig("fig/field_ring.png")
    print("saved fig/field_ring.png")
