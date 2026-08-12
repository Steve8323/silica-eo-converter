#!/usr/bin/env python
"""Directional-coupler design table: L_pi(gap) at 1550 for the w=0.991 core.

Exact supermode calculation (transfer matrix).  Also evaluates the 50/50
splitting-ratio flatness across the 26.8 THz pulse band for each gap —
shorter couplers (smaller gap) are broader-band AND smaller: both design
goals point the same way, bounded below by fab gap limits.
"""
import json
import numpy as np
from multilayer import coupling_length
from dispersion import n_silica, C

W = 0.991
LAM0 = 1.55
BAND_THZ = 26.8  # pulse bandwidth to cover

def ratio_vs_freq(gap, Lc_half, nus):
    """cross-coupled fraction sin^2(pi/2 * L/L_pi(nu)) for L=Lc_half."""
    out = []
    for nu in nus:
        lam = C / nu * 1e6
        Lpi = coupling_length(lam, W, gap, n_silica(lam))
        out.append(np.sin(np.pi / 2 * Lc_half / (Lpi / 1.0)) ** 2
                   if Lpi else np.nan)
    return np.array(out)

if __name__ == "__main__":
    nu0 = C / (LAM0 * 1e-6)
    nus = nu0 + np.linspace(-BAND_THZ / 2, BAND_THZ / 2, 21) * 1e12
    rows = []
    print(f"{'gap[um]':>8} {'L_pi[um]':>9} {'L50[um]':>8} {'ratio min-max over band':>24}")
    for gap in (0.20, 0.30, 0.40, 0.50, 0.70, 1.00):
        Lpi = coupling_length(LAM0, W, gap, n_silica(LAM0))
        L50 = Lpi / 2
        r = ratio_vs_freq(gap, L50, nus)
        rows.append(dict(gap=gap, Lpi=Lpi, L50=L50,
                         rmin=float(np.nanmin(r)), rmax=float(np.nanmax(r))))
        print(f"{gap:8.2f} {Lpi:9.1f} {L50:8.1f}        {np.nanmin(r):.3f} - {np.nanmax(r):.3f}")
    with open("out/coupler_design.json", "w") as f:
        json.dump(rows, f, indent=1)
