#!/usr/bin/env python
"""Master-MRR chi3 demo: degenerate four-wave mixing in the silica ring.

Two CW pumps on adjacent ring resonances (m, m+1) are injected through the
bus.  Silica's chi3 mixes them: new lines appear at 2*f1-f2 and 2*f2-f1 —
the first comb sidebands of a Kerr comb.  chi3 is ARTIFICIALLY BOOSTED so
the effect is visible in an FDTD-scale run; the honest scaling back to real
silica (chi3 = 2.5e-22 m^2/V^2, Q = 1e7-1e9 wedge resonators, mW-class
parametric threshold  P_th ~ 1.54 (pi/2) n^2 V / (n2 lam Q^2)) is in the
report.  What FDTD establishes here: the ring geometry supports the process
and the sidebands land on the neighboring resonances.

usage: python mrr_comb.py F1 F2 CHI3 [runtime]
"""
import json
import sys
import numpy as np
import meep as mp

F1 = float(sys.argv[1])
F2 = float(sys.argv[2])
CHI3 = float(sys.argv[3])
RUNTIME = float(sys.argv[4]) if len(sys.argv) > 4 else 3000.0

R, GAP, W = 8.0, 0.40, 0.991
NSIL = 1.444024
RES = 20

pad = 3.0
sx = 2 * (R + W) + 2 * pad + 4
sy = 2 * (R + W) + 2 * pad + GAP + W + 2
ring_cy = -sy / 2 + pad + W + R + W / 2
bus_y = ring_cy + R + W / 2 + GAP + W / 2

sil_nl = mp.Medium(index=NSIL, chi3=CHI3)
geom = [
    mp.Cylinder(radius=R + W / 2, center=mp.Vector3(0, ring_cy), material=sil_nl),
    mp.Cylinder(radius=R - W / 2, center=mp.Vector3(0, ring_cy),
                material=mp.Medium(index=1.0)),
    mp.Block(size=mp.Vector3(mp.inf, W, mp.inf), center=mp.Vector3(0, bus_y),
             material=sil_nl),
]
src_x = -sx / 2 + 1.8
mon_x = sx / 2 - 1.8
sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom,
                    resolution=RES, boundary_layers=[mp.PML(1.5)],
                    default_material=mp.Medium(index=1.0))
AMP = 1.0
sim.sources = [
    mp.EigenModeSource(mp.ContinuousSource(F1, width=100), amplitude=AMP,
                       center=mp.Vector3(src_x, bus_y), size=mp.Vector3(0, 2.4),
                       eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True),
    mp.EigenModeSource(mp.ContinuousSource(F2, width=100), amplitude=AMP,
                       center=mp.Vector3(src_x, bus_y), size=mp.Vector3(0, 2.4),
                       eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True),
]
# time series at the through port -> spectrum with FWM sidebands
ts = []
def rec(sim):
    ts.append(sim.get_field_point(mp.Ez, mp.Vector3(mon_x, bus_y)).real)
dt_rec = 0.5 / F2 / 8  # 8 samples per optical period min
sim.run(mp.at_every(dt_rec, rec), until=RUNTIME)

t = np.arange(len(ts)) * dt_rec
sig = np.array(ts)
# discard transient (first 40%)
n0 = int(0.4 * len(sig))
sig2 = sig[n0:] * np.hanning(len(sig) - n0)
sp = np.abs(np.fft.rfft(sig2)) ** 2
fr = np.fft.rfftfreq(len(sig2), d=dt_rec)
keep = (fr > 0.55) & (fr < 0.75)
with open(f"out/mrr_comb_chi{CHI3:g}.json", "w") as f:
    json.dump(dict(F1=F1, F2=F2, chi3=CHI3, runtime=RUNTIME,
                   freq=list(fr[keep]), psd=list(sp[keep])), f)
# report line powers at expected comb positions
for name, fx in [("pump1", F1), ("pump2", F2),
                 ("2f1-f2", 2 * F1 - F2), ("2f2-f1", 2 * F2 - F1)]:
    i = np.argmin(np.abs(fr - fx))
    w = 3
    p = sp[max(i - w, 0):i + w + 1].max()
    print(f"{name}: f={fx:.5f}  P={p:.4g}")
