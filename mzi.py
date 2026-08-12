#!/usr/bin/env python
"""2D silica/air MZI with directional couplers only (no Y junctions).

Corrected geometry: ports enter at wide pitch D_A=6 um and s-bend into each
coupler, so nothing couples at the boundaries and the eigenmode launch is
clean.  Centerline of the top guide (bottom is the mirror image):

 [LIN @D_A/2][LB down to D_C/2][LC @D_C/2][LB up][LP phase @D_A/2]
 [LB down][LC][LB up][LIN]

The s-bend approaches add evanescent coupling on top of the straight coupler
section, so the effective 50/50 length is CALIBRATED in FDTD: mode
'coupler' simulates in->DC1->separated arms only and reports the splitting.

Signed output  s = (P_bar - P_cross)/(P_bar + P_cross)  in [-1, 1].

usage: python mzi.py DN LC [coupler|broadband]
"""
import json
import sys
import numpy as np
import meep as mp

DN = float(sys.argv[1])
LC = float(sys.argv[2])
MODE = sys.argv[3] if len(sys.argv) > 3 else "narrow"

W = 0.991
NSIL = 1.444024
RES = 25
GAP_C = 0.30
LB = 14.0
LP = 24.0
D_C = W + GAP_C
D_A = 6.0
PADX = 4.0
LIN = 3.0

COUPLER_ONLY = MODE == "coupler"
seg = [("in", LIN), ("b1", LB), ("c1", LC), ("b2", LB)]
if COUPLER_ONLY:
    seg += [("out", 8.0)]
else:
    seg += [("ph", LP), ("b3", LB), ("c2", LC), ("b4", LB), ("out", LIN)]
edges = np.cumsum([0] + [s[1] for s in seg])
LTOT = edges[-1]
segx = {name: (edges[i], edges[i + 1]) for i, (name, _) in enumerate(seg)}

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
            return D_C / 2 + (D_A - D_C) / 2 * c  # b2, b4
    return D_A / 2

def material(p):
    x = min(max(p.x + LTOT / 2, 0.0), LTOT)
    yc = y_top(x)
    if abs(p.y - yc) <= W / 2:
        dn = DN if (not COUPLER_ONLY and segx["ph"][0] <= x <= segx["ph"][1]) else 0.0
        return mp.Medium(index=NSIL + dn)
    if abs(p.y + yc) <= W / 2:
        return mp.Medium(index=NSIL)
    return mp.Medium(index=1.0)

fcen = 1 / 1.55
df = (0.35 if MODE == "broadband" else 0.10) * fcen
sx = LTOT + 2 * PADX
sy = D_A + W + 2 * 4.0
cell = mp.Vector3(sx, sy)
src_x = -sx / 2 + PADX * 0.55
mon_x = sx / 2 - PADX * 0.55

sim = mp.Simulation(cell_size=cell, resolution=RES,
                    boundary_layers=[mp.PML(1.5)],
                    material_function=material,
                    extra_materials=[mp.Medium(index=NSIL + DN + 1e-4)],
                    default_material=mp.Medium(index=1.0))
sim.sources = [mp.EigenModeSource(
    mp.GaussianSource(fcen, fwidth=df),
    center=mp.Vector3(src_x, D_A / 2), size=mp.Vector3(0, 3.6),
    eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]

nfreq = 161
fl_bar = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(
    center=mp.Vector3(mon_x, D_A / 2), size=mp.Vector3(0, 3.6)))
fl_cross = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(
    center=mp.Vector3(mon_x, -D_A / 2), size=mp.Vector3(0, 3.6)))
fl_in = sim.add_flux(fcen, df, nfreq, mp.FluxRegion(
    center=mp.Vector3(src_x + 1.2, D_A / 2), size=mp.Vector3(0, 3.6)))

sim.run(until_after_sources=mp.stop_when_fields_decayed(
    60, mp.Ez, mp.Vector3(mon_x, D_A / 2), 1e-8))

freqs = np.array(mp.get_flux_freqs(fl_bar))
Pb = np.array(mp.get_fluxes(fl_bar))
Pc = np.array(mp.get_fluxes(fl_cross))
Pi = np.array(mp.get_fluxes(fl_in))
s = (Pb - Pc) / (Pb + Pc)
tag = f"{DN:g}_LC{LC:g}_{MODE}"
with open(f"out/mzi_{tag}.json", "w") as f:
    json.dump(dict(dn=DN, LC=LC, mode=MODE, lam_nm=list(1000 / freqs),
                   P_bar=list(Pb), P_cross=list(Pc), P_in=list(Pi),
                   s=list(s)), f)
i0 = np.argmin(np.abs(freqs - fcen))
if COUPLER_ONLY:
    print(f"LC={LC:g}: cross fraction={Pc[i0]/(Pb[i0]+Pc[i0]):.4f} "
          f"(bar={Pb[i0]:.4g}, cross={Pc[i0]:.4g}, in={Pi[i0]:.4g})")
else:
    print(f"DN={DN:g} LC={LC:g}: s(1550)={s[i0]:+.4f} bar={Pb[i0]:.4g} "
          f"cross={Pc[i0]:.4g} sum/in={(Pb[i0]+Pc[i0])/Pi[i0]:.3f}")
