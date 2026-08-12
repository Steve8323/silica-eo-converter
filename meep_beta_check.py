#!/usr/bin/env python
"""MEEP FDTD validation of the analytic slab dispersion.

2D silica strip (dispersive fused quartz from the MEEP materials library) in
air.  Launches a broadband TE (Ez) eigenmode pulse, records complex DFT
fields at two cross-sections dx apart, and extracts

    ng(nu)    = d(dphi)/dom / dx
    beta2(nu) = d^2(dphi)/dom^2 / dx

Derivatives of the phase DIFFERENCE are free of the 2*pi*m branch ambiguity,
so this is a genuine, prior-free check of the analytic ng and beta2.
"""
import json
import numpy as np
import meep as mp
from meep.materials import fused_quartz

W = 0.991           # core width (um) — zero-GVD design
DX = 200.0          # monitor separation (um)
RES = 25
CONST = __import__("sys").argv[1:2] == ["const"]
LAM0 = 1.55
DF_FRAC = 0.35      # fractional bandwidth of source

def run():
    fcen = 1 / LAM0
    df = DF_FRAC * fcen
    pad = 20.0
    sx = DX + 2 * pad
    sy = 8.0
    cell = mp.Vector3(sx, sy)
    mat = mp.Medium(index=1.444024) if CONST else fused_quartz
    geom = [mp.Block(size=mp.Vector3(mp.inf, W, mp.inf), center=mp.Vector3(),
                     material=mat)]
    src = [mp.EigenModeSource(mp.GaussianSource(fcen, fwidth=df),
                              center=mp.Vector3(-sx / 2 + 8, 0),
                              size=mp.Vector3(0, sy - 2),
                              eig_band=1, eig_parity=mp.ODD_Z + mp.EVEN_Y,
                              eig_match_freq=True)]
    sim = mp.Simulation(cell_size=cell, geometry=geom, sources=src,
                        resolution=RES, boundary_layers=[mp.PML(1.5)],
                        symmetries=[mp.Mirror(mp.Y, phase=1)],
                        Courant=0.5 if CONST else 0.3,
                        default_material=mp.Medium(index=1.0))
    nfreq = 401
    x1, x2 = -DX / 2, DX / 2
    mon = [sim.add_dft_fields([mp.Ez], fcen, df * 1.2, nfreq,
                              center=mp.Vector3(x, 0),
                              size=mp.Vector3(0, 0))
           for x in (x1, x2)]
    # pulse needs ~ ng*(DX+2*pad) ~ 360 meep-time to clear the far monitor.
    # Hard cap: the guided pulse passes both monitors well before t=500; the
    # residual is non-guided glancing radiation that decays for thousands of
    # time units without touching the guided-mode DFT phases.
    sim.run(until_after_sources=500)

    f1 = np.array([sim.get_dft_array(mon[0], mp.Ez, i) for i in range(nfreq)]).ravel()
    f2 = np.array([sim.get_dft_array(mon[1], mp.Ez, i) for i in range(nfreq)]).ravel()
    freqs = np.linspace(fcen - 0.6 * df, fcen + 0.6 * df, nfreq)
    dphi = np.unwrap(np.angle(f2 / f1))
    keep = (np.abs(f1) > 1e-3 * np.abs(f1).max())
    return freqs, dphi, keep, np.abs(f1)

if __name__ == "__main__":
    freqs, dphi, keep, amp = run()
    # meep units: length um, freq c/um. Convert to SI.
    C = 299792458.0
    nus = freqs * C / 1e-6            # Hz
    om = 2 * np.pi * nus
    dx_m = DX * 1e-6
    # smooth polynomial fit of dphi(om) over the kept band
    o, p_ = om[keep], dphi[keep]
    pc = np.polyfit(o - o.mean(), p_, 7)
    d1 = np.polyder(pc, 1)
    d2 = np.polyder(pc, 2)
    ng = np.polyval(d1, o - o.mean()) / dx_m * C
    b2 = np.polyval(d2, o - o.mean()) / dx_m
    out = dict(nu=list(nus[keep]), ng=list(ng), beta2=list(b2),
               amp=list(amp[keep] / amp.max()), const=CONST)
    with open("out/meep_beta_check.json", "w") as f:
        json.dump(out, f)
    i0 = np.argmin(np.abs(nus[keep] - C / 1.55e-6))
    print(f"MEEP  @1550: ng={ng[i0]:.4f}, beta2={b2[i0]*1e27:.2f} ps^2/km")
    print(f"band: {nus[keep].min()/1e12:.1f}-{nus[keep].max()/1e12:.1f} THz, "
          f"{keep.sum()} pts")
