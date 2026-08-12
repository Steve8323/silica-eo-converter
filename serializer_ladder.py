#!/usr/bin/env python
"""Counter-directional tap-ladder serializer (FIR -> exactly blank tail).

Bus A (top, flows right) carries the input pulse.  N connector arms tap it
with ladder ratios s_k = 1/(N+1-k) (equal 1/N energy per arm), U-turn, and
merge onto bus B (bottom, flows LEFT) with ratios m_k = 1/k, giving uniform
slot transmission T_k = 1/N^2 — the static-passive optimum for crisp slots.
Slot delay of tap k = 2*x_k*ng/c: set by tap POSITION, so uneven spacing is
free (pitch list below is deliberately non-uniform).  After the merge the
residual in each arm bends 90 deg down, crosses bus B perpendicularly
(clean, <~1% scatter) and dumps into the PML.

No waveguide merges anywhere — every junction is an evanescent coupler.

usage: python serializer_ladder.py N [res]
"""
import json
import sys
import numpy as np
import meep as mp

N = int(sys.argv[1])
RES = int(sys.argv[2]) if len(sys.argv) > 2 else 15

W = 0.991
NSIL = 1.444024
GAP = 0.12
LPI = 11.31           # full-transfer length at gap 0.12 (supermode)
RHO = 6.0             # U-bend radius: 3.5 radiated ~10%/half-turn (ring Q data); 6.0 -> ~1%
FCEN = 1 / 1.55
FW = 0.15 * FCEN

GAP_M = 0.40    # merge-side gap: lowers the arc-coupling floor below 1/N

# FDTD-calibrated tap curve at gap 0.12 with rho=6 (calib_cell sweeps):
CAL_L = [0.15, 2.5, 5.0, 8.0]
CAL_K = [0.434, 0.713, 0.949, 1.253]      # kappa(L1) [rad], monotone

def L_tap(ratio):
    k = np.arcsin(np.sqrt(min(ratio, 1.0)))
    if k <= CAL_K[0]:
        return 0.15
    if k > CAL_K[-1]:                       # linear extrapolation at top
        return CAL_L[-1] + (k - CAL_K[-1]) / 0.101
    return float(np.interp(k, CAL_K, CAL_L))

def L_mrg(ratio):
    # extrapolated merge curve at gap 0.40: kappa = 0.45 + 0.075*L
    k = np.arcsin(np.sqrt(min(ratio, 1.0)))
    return max((k - 0.45) / 0.075, 0.15)

PITCH = [16, 18, 20, 24, 26, 28, 30, 32][: N]      # uneven, front-loaded
# split: s_k = 1/(N-k) leaves 1/N per arm; merge MIRRORS it (full transfer
# at the LAST tap, which no other slot passes): T_k = 1/N^2 uniform.
s_ratio = [1.0 / (N - k) for k in range(N)]
m_ratio = [1.0 / (N - k) for k in range(N)]
L1 = [L_tap(r) for r in s_ratio]
L2 = [L_mrg(r) for r in m_ratio]
if N == 6:
    # trim iteration 1: flatten measured slot ripple
    # (slots were 0.48,0.71,0.72,1.00,0.32,0.22 -> cut hot mid taps/merges,
    #  feed the starved late slots via the (1-s) chain)
    # trim iteration 2 (iter-1 slots: 1.00,0.34,0.36,0.30,0.39,0.21 —
    #  m_1 boost overshot 3x: merge curve steeper than extrapolated near floor)
    L1 = [0.15, 0.15, 0.50, 0.65, 3.60, 10.5]
    L2 = [0.15, 0.15, 0.15, 0.15, 3.50, 16.0]

xs = np.concatenate([[10.0], 10.0 + np.cumsum(PITCH[:-1])]) if N > 1 else np.array([10.0])
SPAN = xs[-1] + L1[-1] + 2 * RHO + 8
sx = SPAN + 14
# vertical levels
yA = 6.0
y_arm = yA - (W + GAP)          # tap-arm level (below bus A)
y_mrg = y_arm - 2 * RHO         # merge-arm level
yB = y_mrg - (W + GAP_M)        # bus B
sy = (yA + 5) - (yB - 7)
y0 = (yA + 5 + yB - 7) / 2      # cell center offset

segs = []   # ((x1,y1),(x2,y2)) straight centerlines
arcs = []   # (cx, cy, r, th1, th2) centerline arcs
segs.append(((-1e3, yA), (1e3, yA)))            # bus A
segs.append(((-1e3, yB), (1e3, yB)))            # bus B
for k in range(N):
    x = xs[k]
    xu = x + L1[k]              # U-bend start
    segs.append(((x - 1.0, y_arm), (xu, y_arm)))            # tap section
    arcs.append((xu, y_arm - RHO, RHO, -90, 90))            # right semicircle
    xm = xu
    segs.append(((xm - L2[k], y_mrg), (xm, y_mrg)))         # merge section
    xd = xm - L2[k] - 3.0                                   # down-bend corner
    segs.append(((xd, y_mrg), (xm - L2[k], y_mrg)))         # short lead
    arcs.append((xd, y_mrg - 4.0, 4.0, 90, 180))            # quarter bend down
    segs.append(((xd - 4.0, y_mrg - 4.0), (xd - 4.0, yB - 8)))  # dump stub

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
            if t1 <= th <= t2 or t1 <= th + 360 <= t2:
                return mp.Medium(index=NSIL)
    return mp.Medium(index=1.0)

sim = mp.Simulation(cell_size=mp.Vector3(sx, sy), resolution=RES,
                    boundary_layers=[mp.PML(1.5)],
                    material_function=material,
                    default_material=mp.Medium(index=1.0))
sim.sources = [mp.EigenModeSource(
    mp.GaussianSource(FCEN, fwidth=FW),
    center=mp.Vector3(-sx / 2 + 2.2, yA - y0), size=mp.Vector3(0, 3.0),
    eig_band=1, eig_parity=mp.ODD_Z, eig_match_freq=True,
    direction=mp.X)]

nf = 81
fl_in = sim.add_flux(FCEN, FW * 1.4, nf, mp.FluxRegion(
    center=mp.Vector3(-sx / 2 + 3.4, yA - y0), size=mp.Vector3(0, 3.0)))
out_p = mp.Vector3(-sx / 2 + 2.2, yB - y0)
fl_out = sim.add_flux(FCEN, FW * 1.4, nf, mp.FluxRegion(
    center=out_p, size=mp.Vector3(0, 3.0)))

ts = []
DT = 0.12
sim.run(mp.at_every(DT, lambda s: ts.append(
    s.get_field_point(mp.Ez, out_p).real)),
    until=1.6 * 2 * SPAN * 1.5 + 260)

sig = np.array(ts)
S = np.fft.fft(sig); S[len(S) // 2:] = 0
env = (2 * np.abs(np.fft.ifft(S))) ** 2
t_fs = np.arange(len(sig)) * DT / 0.299792458
Pin = np.array(mp.get_fluxes(fl_in)).sum()
Pout = -np.array(mp.get_fluxes(fl_out)).sum()   # leftward flux
# +33 um internal arm path (U-bend + leads), calibrated from run 1
slot_fs = [(2 * (x + 4) + 33) * 1.4948 / 0.299792458 for x in xs]
with open(f"out/ladder_N{N}.json", "w") as f:
    json.dump(dict(N=N, xs=list(xs), pitch=PITCH, s=s_ratio, m=m_ratio,
                   L1=L1, L2=L2, slot_fs=slot_fs,
                   t_fs=list(t_fs[::3]), env=list(env[::3]),
                   ratio=float(Pout / Pin)), f)
print(f"ladder N={N}: E_out/E_in={Pout/Pin:.4f} (ideal 1/N={1/N:.3f}); "
      f"slots(fs)={[f'{s:.0f}' for s in slot_fs]}")
