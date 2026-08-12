#!/usr/bin/env python
"""Calibrate the ladder's tap and merge couplers IN SITU (with their arcs).

Single connector cell: bus A + arm (tap section L1, U-bend, merge section L2,
dump stub) + bus B.  Measures
    s = 1 - P_A_through / P_in          (tap ratio, incl. arc contribution)
    m = P_B_out / (P_in - P_A_through)  (merge ratio of what reached the arm)

usage: python calib_cell.py L1 L2
"""
import json
import sys
import numpy as np
import meep as mp

L1 = float(sys.argv[1])
L2 = float(sys.argv[2])
W, NSIL, GAP, RES = 0.991, 1.444024, 0.12, 15
GAP_M = float(sys.argv[3]) if len(sys.argv) > 3 else 0.30
RHO = 6.0
FCEN, FW = 1 / 1.55, 0.15 / 1.55

x_tap = 10.0
xu = x_tap + L1
SPAN = xu + 2 * RHO + 8
sx = SPAN + 14
yA = 6.0
y_arm = yA - (W + GAP)
y_mrg = y_arm - 2 * RHO
yB = y_mrg - (W + GAP_M)
sy = (yA + 5) - (yB - 7)
y0 = (yA + 5 + yB - 7) / 2

segs = [((-1e3, yA), (1e3, yA)), ((-1e3, yB), (1e3, yB)),
        ((x_tap - 1.0, y_arm), (xu, y_arm)),
        ((xu - L2, y_mrg), (xu, y_mrg))]
xd = xu - L2 - 2.0
segs += [((xd, y_mrg), (xu - L2, y_mrg)),
         ((xd - 2.0, y_mrg - 2.0), (xd - 2.0, yB - 8))]
arcs = [(xu, y_arm - RHO, RHO, -90, 90), (xd, y_mrg - 2.0, 2.0, 90, 180)]

def dist_seg(p, a, b):
    ax, ay = a; bx, by = b
    t = ((p[0] - ax) * (bx - ax) + (p[1] - ay) * (by - ay)) / \
        max((bx - ax) ** 2 + (by - ay) ** 2, 1e-12)
    t = min(max(t, 0.0), 1.0)
    return np.hypot(p[0] - (ax + t * (bx - ax)), p[1] - (ay + t * (by - ay)))

def material(p):
    x, y = p.x + sx / 2, p.y + y0
    for a, b in segs:
        if dist_seg((x, y), a, b) <= W / 2:
            return mp.Medium(index=NSIL)
    for cx, cy, r, t1, t2 in arcs:
        d = np.hypot(x - cx, y - cy)
        if abs(d - r) <= W / 2:
            th = np.degrees(np.arctan2(y - cy, x - cx))
            if t1 <= th <= t2:
                return mp.Medium(index=NSIL)
    return mp.Medium(index=1.0)

sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), resolution=RES,
                    boundary_layers=[mp.PML(1.5)],
                    material_function=material,
                    default_material=mp.Medium(index=1.0))
sim.sources = [mp.EigenModeSource(
    mp.GaussianSource(FCEN, fwidth=FW),
    center=mp.Vector3(-sx / 2 + 2.2, yA - y0), size=mp.Vector3(0, 3.0),
    eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]
nf = 41
def flux(x, y):
    return sim.add_flux(FCEN, FW * 1.2, nf, mp.FluxRegion(
        center=mp.Vector3(x, y), size=mp.Vector3(0, 3.0)))
fl_in = flux(-sx / 2 + 3.4, yA - y0)
fl_thr = flux(sx / 2 - 2.5, yA - y0)
fl_out = flux(-sx / 2 + 2.2, yB - y0)
sim.run(until=2.2 * SPAN * 1.5 + 260)
Pin = np.array(mp.get_fluxes(fl_in)).sum()
Pthr = np.array(mp.get_fluxes(fl_thr)).sum()
Pout = -np.array(mp.get_fluxes(fl_out)).sum()
s = 1 - Pthr / Pin
m = Pout / max(Pin - Pthr, 1e-12)
with open(f"out/calib_L1_{L1:g}_L2_{L2:g}.json", "w") as f:
    json.dump(dict(L1=L1, L2=L2, s=float(s), m=float(m)), f)
print(f"L1={L1:g} L2={L2:g}: s={s:.4f} m={m:.4f}")
