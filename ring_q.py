#!/usr/bin/env python
"""Intrinsic (bend-loss) Q of a silica/air 2D ring, one radius per run.

usage: python ring_q.py R_um

Ring = annulus of silica (n=1.444024, constant — narrowband sim) with the
w=0.991 um zero-GVD cross-section.  Point Ez dipole inside the guide excites
whispering-gallery modes near 1550 nm; Harminv extracts (f, Q).  This Q is
radiation-limited only (no material/roughness loss in 2D FDTD), i.e. the
bend-loss ceiling that sets the minimum ring/bend radius of the circuit.
"""
import json
import sys
import meep as mp

R = float(sys.argv[1])
DN_ARG = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0
W = 0.991
N_SILICA = 1.444024
RES = 30
FCEN = 1 / 1.55
DF = 0.08 * FCEN

pad = 3.0
s = 2 * (R + W + pad)
cell = mp.Vector3(s, s)
geom = [mp.Cylinder(radius=R + W / 2, material=mp.Medium(index=N_SILICA + DN_ARG)),
        mp.Cylinder(radius=R - W / 2, material=mp.Medium(index=1.0))]
src = [mp.Source(mp.GaussianSource(FCEN, fwidth=DF), component=mp.Ez,
                 center=mp.Vector3(R, 0.11))]
sim = mp.Simulation(cell_size=cell, geometry=geom, sources=src,
                    resolution=RES, boundary_layers=[mp.PML(1.5)])
h = mp.Harminv(mp.Ez, mp.Vector3(-R, 0.05), FCEN, DF)
sim.run(mp.after_sources(h), until_after_sources=1200)

rows = [dict(f=m.freq, lam_nm=1000 / m.freq, Q=m.Q, amp=abs(m.amp))
        for m in h.modes if m.Q > 0]
rows.sort(key=lambda r: -r["amp"])
with open(f"out/ring_q_R{R:g}_dn{DN_ARG:g}.json", "w") as f:
    json.dump(dict(R=R, modes=rows), f, indent=1)
for r in rows[:6]:
    print(f"R={R:g}: lam={r['lam_nm']:.2f} nm  Q={r['Q']:.3g}  amp={r['amp']:.2g}")
