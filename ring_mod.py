#!/usr/bin/env python
"""Bus-coupled silica MRR: transmission spectrum vs electrode Delta-n.

usage: python ring_mod.py R GAP DN [res]

The ring index is NSIL + DN (electrode surrogate: uniform coverage of the
ring by the field of coplanar electrodes; the real Delta-n <-> V map is in
eo_mapping.py).  Broadband source in the bus; DFT flux at the through port
resolves the resonance dips.  Runs until the through-port fields decay, so
high-Q dips are fully resolved.
"""
import json
import sys
import numpy as np
import meep as mp

R = float(sys.argv[1])
GAP = float(sys.argv[2])
DN = float(sys.argv[3])
RES = int(sys.argv[4]) if len(sys.argv) > 4 else 20

W = 0.991
NSIL = 1.444024
fcen = 1 / 1.55
df = 0.06 * fcen

pad = 3.0
sx = 2 * (R + W) + 2 * pad + 4
sy = 2 * (R + W) + 2 * pad + GAP + W + 2
ring_cy = -sy / 2 + pad + W + R + W / 2          # ring center y
bus_y = ring_cy + R + W / 2 + GAP + W / 2        # bus centerline

geom = [
    mp.Cylinder(radius=R + W / 2, center=mp.Vector3(0, ring_cy),
                material=mp.Medium(index=NSIL + DN)),
    mp.Cylinder(radius=R - W / 2, center=mp.Vector3(0, ring_cy),
                material=mp.Medium(index=1.0)),
    mp.Block(size=mp.Vector3(mp.inf, W, mp.inf), center=mp.Vector3(0, bus_y),
             material=mp.Medium(index=NSIL)),
]
src_x = -sx / 2 + 1.8
mon_x = sx / 2 - 1.8
sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom,
                    resolution=RES, boundary_layers=[mp.PML(1.5)],
                    default_material=mp.Medium(index=1.0))
sim.sources = [mp.EigenModeSource(
    mp.GaussianSource(fcen, fwidth=df),
    center=mp.Vector3(src_x, bus_y), size=mp.Vector3(0, 2.4),
    eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]
nfreq = 2001
fl = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(
    center=mp.Vector3(mon_x, bus_y), size=mp.Vector3(0, 2.4)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(
    100, mp.Ez, mp.Vector3(mon_x, bus_y), 1e-7))

freqs = np.array(mp.get_flux_freqs(fl))
P = np.array(mp.get_fluxes(fl))
with open(f"out/ring_mod_R{R:g}_g{GAP:g}_dn{DN:g}.json", "w") as f:
    json.dump(dict(R=R, gap=GAP, dn=DN, res=RES,
                   lam_nm=list(1000 / freqs), P=list(P)), f)
print(f"R={R} gap={GAP} dn={DN:g}: done, t_end={sim.meep_time():.0f}")
