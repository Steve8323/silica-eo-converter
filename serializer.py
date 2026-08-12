#!/usr/bin/env python
"""All-pass ring-chain serializer: N values -> uneven pulse cluster, zero loss.

N rings (different radii -> different round-trip delays, deliberately
non-commensurate) side-coupled to ONE bus.  Each ring is all-pass: every
photon returns to the same bus mode, so adding rings adds ZERO scaling loss
(only 2D bend radiation, ~0.1%/round-trip at R>=6).  Impulse response =
input pulse convolved with each ring's echo train: a quick UNEVENLY-SPACED
pulse cluster at the head of the frame, geometric decay -> blank tail.
Each ring is an electrode site (EO Delta-n shifts its echo phase).

Compressed for a 10-min budget: res 15, short cell, single broadband pulse,
time-domain point record at the output + broadband flux in/out for the
energy-conservation (zero-loss) check.

usage: python serializer.py N [gap]
"""
import json
import sys
import numpy as np
import meep as mp

N = int(sys.argv[1])
GAP = float(sys.argv[2]) if len(sys.argv) > 2 else 0.15
RADII = [6.0, 7.3, 8.1, 9.2][:N]
W = 0.991
NSIL = 1.444024
RES = 15
FCEN = 1 / 1.55
FW = 0.15 * FCEN

pad = 6.0
xs, x = [], 0.0
for R in RADII:
    x += R + 2.5
    xs.append(x)
    x += R + 2.5
LTOT = x + 2 * pad
sx = LTOT
sy = 2 * max(RADII) + W + GAP + 9

bus_y = -sy / 2 + 3.5
geom = [mp.Block(size=mp.Vector3(mp.inf, W, mp.inf), center=mp.Vector3(0, bus_y),
                 material=mp.Medium(index=NSIL))]
for R, xc in zip(RADII, xs):
    cy = bus_y + W / 2 + GAP + R + W / 2
    c = mp.Vector3(xc - sx / 2 + pad, cy)
    geom += [mp.Cylinder(radius=R + W / 2, center=c, material=mp.Medium(index=NSIL)),
             mp.Cylinder(radius=R - W / 2, center=c, material=mp.Medium(index=1.0))]

src_x = -sx / 2 + 2.0
out_x = sx / 2 - 2.0
sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), geometry=geom,
                    resolution=RES, boundary_layers=[mp.PML(1.5)],
                    default_material=mp.Medium(index=1.0))
sim.sources = [mp.EigenModeSource(
    mp.GaussianSource(FCEN, fwidth=FW),
    center=mp.Vector3(src_x, bus_y), size=mp.Vector3(0, 3.0),
    eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True)]

nf = 101
fl_in = sim.add_flux(FCEN, FW * 1.5, nf, mp.FluxRegion(
    center=mp.Vector3(src_x + 1.0, bus_y), size=mp.Vector3(0, 3.0)))
fl_out = sim.add_flux(FCEN, FW * 1.5, nf, mp.FluxRegion(
    center=mp.Vector3(out_x, bus_y), size=mp.Vector3(0, 3.0)))

ts = []
DT = 0.12
sim.run(mp.at_every(DT, lambda s: ts.append(
    s.get_field_point(mp.Ez, mp.Vector3(out_x, bus_y)).real)),
    until=260 + 120 * N + 500)  # +500: capture the full geometric echo tail

Pin = np.array(mp.get_fluxes(fl_in)).sum()
Pout = np.array(mp.get_fluxes(fl_out)).sum()
sig = np.array(ts)
# envelope: |analytic signal| via FFT (cheap Hilbert)
S = np.fft.fft(sig)
S[len(S) // 2:] = 0
env = 2 * np.abs(np.fft.ifft(S))
t_fs = np.arange(len(sig)) * DT * 1e-6 / 299792458.0 * 1e15  # meep t -> fs

rts = [2 * np.pi * R * 1.4948 / 0.299792458 for R in RADII]  # fs round trips
with open(f"out/serializer_N{N}_g{GAP:g}.json", "w") as f:
    json.dump(dict(N=N, gap=GAP, radii=RADII, rt_fs=rts,
                   t_fs=list(t_fs[::3]), env=list((env ** 2)[::3]),
                   Ein=float(Pin), Eout=float(Pout),
                   ratio=float(Pout / Pin)), f)
print(f"N={N} gap={GAP}: E_out/E_in = {Pout/Pin:.4f}  "
      f"round-trips(fs) = {[f'{r:.0f}' for r in rts]}")
