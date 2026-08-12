#!/usr/bin/env python
"""Exact linear propagation over L=5 mm using the full slab beta(omega).

No Taylor expansion: E_out(om) = E_in(om) * exp(i beta(om) L) with beta from
the exact dispersion solver (Malitson silica, air clad).  Binary-searches the
shortest Gaussian T0 whose FWHM broadening factor stays <= sqrt(2).

Widths tested: w = 0.744 um (single-mode limit) and w = 0.991 um (zero-GVD).
"""
import json
import numpy as np
from dispersion import slab_beta_TE0, C

L = 5e-3
LAM0 = 1.55e-6
NU0 = C / LAM0

def beta_grid(w_um, span_thz=120.0, n=4001):
    nus = NU0 + np.linspace(-span_thz / 2, span_thz / 2, n) * 1e12
    lams = C / nus * 1e6
    b = np.array([slab_beta_TE0(l, w_um) for l in lams])
    return nus, b

def fwhm(t, I):
    I = I / I.max()
    above = np.where(I >= 0.5)[0]
    if len(above) < 2:
        return np.nan
    i0, i1 = above[0], above[-1]
    # linear interp at the edges
    def cross(i, j):
        return t[i] + (0.5 - I[i]) * (t[j] - t[i]) / (I[j] - I[i])
    t_lo = cross(i0 - 1, i0) if i0 > 0 else t[i0]
    t_hi = cross(i1 + 1, i1) if i1 < len(t) - 1 else t[i1]
    return t_hi - t_lo

def broadening(T0, nus, beta, mask):
    """FWHM_out / FWHM_in for Gaussian T0 through exp(i beta L)."""
    om = 2 * np.pi * nus
    om0 = 2 * np.pi * NU0
    # Gaussian E(om) ~ exp(-T0^2 (om-om0)^2 / 2)
    Ein = np.exp(-T0 ** 2 * (om - om0) ** 2 / 2)
    if Ein[mask][0] > 1e-3 or Ein[mask][-1] > 1e-3:
        return None  # spectrum falls outside guided/valid band
    # remove linear phase (group delay) for a centered time window
    p = np.polyfit(om - om0, beta, 1)
    phi = (beta - np.polyval(p, om - om0)) * L
    Eout = Ein * np.exp(1j * phi)
    # time grid via FFT
    dom = om[1] - om[0]
    t = np.fft.fftshift(np.fft.fftfreq(len(om), d=dom / (2 * np.pi)))
    Et_in = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Ein)))
    Et_out = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Eout)))
    return fwhm(t, np.abs(Et_out) ** 2) / fwhm(t, np.abs(Et_in) ** 2)

def shortest(nus, beta):
    mask = ~np.isnan(beta)
    # restrict everything to valid band
    nus2, beta2 = nus[mask], beta[mask]
    m2 = np.ones(len(nus2), bool)
    lo, hi = 3e-15, 300e-15  # T0 search range
    for _ in range(60):
        mid = np.sqrt(lo * hi)
        b = broadening(mid, nus2, beta2, m2)
        if b is None or b > np.sqrt(2):
            lo = mid
        else:
            hi = mid
    return hi

if __name__ == "__main__":
    res = {}
    for w in (0.744, 0.991):
        nus, b = beta_grid(w)
        T0 = shortest(nus, b)
        fw = 1.665 * T0
        bw = 0.44 / fw
        res[str(w)] = dict(T0=T0, fwhm=fw, bw_thz=bw / 1e12)
        print(f"w={w} um: shortest T0={T0*1e15:.1f} fs, FWHM={fw*1e15:.1f} fs, "
              f"transform-limited BW={bw/1e12:.1f} THz  (exact beta(omega), L=5mm)")
    with open("out/pulse_5mm.json", "w") as f:
        json.dump(res, f, indent=1)
