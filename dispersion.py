#!/usr/bin/env python
"""Silica-slab (2D) waveguide dispersion — EXACT.

Symmetric slab: silica core (Malitson Sellmeier), air cladding, TE
polarization (E out of plane = meep Ez).  Solves the transcendental slab
equation per frequency, builds beta(omega), extracts beta1/beta2/beta3 by
finite differences, then finds

  * the core width w* that minimizes |beta2| at 1550 nm (zero-GVD point)
  * the shortest Gaussian pulse T0 that keeps broadening <= sqrt(2) over
    L = 5 mm, including the beta3 limit:
        beta2 limit : T0 >= sqrt(|beta2| L)
        beta3 limit : T0 >= (|beta3| L)^(1/3)

No approximations beyond the 2D (slab) reduction itself, which matches the
2D MEEP model used everywhere else in this project.
"""
import json
import numpy as np
from scipy.optimize import brentq

C = 299792458.0  # m/s

def n_silica(lam_um):
    """Malitson fused-silica Sellmeier (lam in um)."""
    l2 = lam_um ** 2
    n2 = (1
          + 0.6961663 * l2 / (l2 - 0.0684043 ** 2)
          + 0.4079426 * l2 / (l2 - 0.1162414 ** 2)
          + 0.8974794 * l2 / (l2 - 9.896161 ** 2))
    return np.sqrt(n2)

def slab_beta_TE0(lam_um, w_um, n_clad=1.0, n1=None):
    """Exact propagation constant beta [rad/m] of the fundamental even TE
    slab mode: tan(kappa d/2) = gamma/kappa."""
    if n1 is None:
        n1 = n_silica(lam_um)
    k0 = 2 * np.pi / (lam_um * 1e-6)
    d = w_um * 1e-6

    def f(neff):
        kap = k0 * np.sqrt(n1 ** 2 - neff ** 2)
        gam = k0 * np.sqrt(neff ** 2 - n_clad ** 2)
        return np.tan(kap * d / 2) - gam / kap

    lo, hi = n_clad + 1e-9, n1 - 1e-9
    # TE0 root is the first branch: kappa*d/2 in (0, pi/2). Restrict neff so
    # kappa*d/2 < pi/2 to stay on that branch.
    kap_max = np.pi / d  # kappa*d/2 = pi/2
    neff_min_branch = np.sqrt(max(n1 ** 2 - (kap_max / k0) ** 2, n_clad ** 2))
    lo = max(lo, neff_min_branch + 1e-12)
    if f(lo) * f(hi) > 0:  # mode cut off / numeric edge
        return np.nan
    neff = brentq(f, lo, hi, xtol=1e-14)
    return neff * k0

def beta_derivs(w_um, lam0_um=1.55, dnu=5e9, span=201):
    """beta1..beta3 at lam0 by central finite differences on a nu grid."""
    nu0 = C / (lam0_um * 1e-6)
    nus = nu0 + dnu * np.arange(-(span // 2), span // 2 + 1)
    lams = C / nus * 1e6
    betas = np.array([slab_beta_TE0(l, w_um) for l in lams])
    if np.any(np.isnan(betas)):
        return None
    om = 2 * np.pi * nus
    # fit local polynomial (degree 5) for smooth derivatives
    p = np.polyfit(om - 2 * np.pi * nu0, betas, 5)
    b1 = p[-2]
    b2 = 2 * p[-3]
    b3 = 6 * p[-4]
    return dict(neff=betas[span // 2] / (2 * np.pi * nu0 / C),
                ng=b1 * C, beta2=b2, beta3=b3)

def shortest_pulse(beta2, beta3, L=5e-3):
    """Shortest Gaussian T0 (1/e half-width, s) with broadening<=sqrt(2) over L."""
    t2 = np.sqrt(abs(beta2) * L)
    t3 = abs(beta3 * L) ** (1 / 3.0)
    return max(t2, t3), t2, t3

if __name__ == "__main__":
    ws = np.linspace(0.30, 3.0, 111)
    rows = []
    for w in ws:
        d = beta_derivs(w)
        if d is None:
            continue
        T0, t2, t3 = shortest_pulse(d["beta2"], d["beta3"])
        rows.append(dict(w=w, **d, T0=T0, T0_b2=t2, T0_b3=t3))
    # zero-GVD width: minimize |beta2|
    best = min(rows, key=lambda r: abs(r["beta2"]))
    # refine around it
    wfine = np.linspace(best["w"] - 0.06, best["w"] + 0.06, 61)
    for w in wfine:
        d = beta_derivs(w)
        if d is None:
            continue
        T0, t2, t3 = shortest_pulse(d["beta2"], d["beta3"])
        rows.append(dict(w=w, **d, T0=T0, T0_b2=t2, T0_b3=t3))
    rows.sort(key=lambda r: r["w"])
    best = min(rows, key=lambda r: abs(r["beta2"]))

    out = dict(rows=rows, best=best,
               bulk_n=n_silica(1.55),
               note="beta2 [s^2/m], beta3 [s^3/m], T0 [s] Gaussian 1/e half-width; "
                    "FWHM = 1.665*T0; L = 5 mm; broadening criterion sqrt(2)")
    with open("out/dispersion.json", "w") as f:
        json.dump(out, f, indent=1)

    print(f"bulk silica n(1550) = {n_silica(1.55):.6f}")
    print(f"{'w[um]':>7} {'neff':>8} {'ng':>7} {'beta2[ps^2/km]':>15} {'beta3[ps^3/km]':>15} {'T0[fs]':>8} {'FWHM[fs]':>9}")
    for r in rows[::10]:
        print(f"{r['w']:7.2f} {r['neff']:8.4f} {r['ng']:7.4f} "
              f"{r['beta2']*1e27:15.2f} {r['beta3']*1e42:15.3f} "
              f"{r['T0']*1e15:8.1f} {r['T0']*1.665e15:9.1f}")
    print("\nZero-GVD (min |beta2|):")
    r = best
    print(f"  w* = {r['w']:.3f} um, neff={r['neff']:.4f}, ng={r['ng']:.4f}, "
          f"beta2={r['beta2']*1e27:.3f} ps^2/km, beta3={r['beta3']*1e42:.3f} ps^3/km")
    print(f"  T0 = {r['T0']*1e15:.1f} fs (b2 limit {r['T0_b2']*1e15:.1f}, b3 limit {r['T0_b3']*1e15:.1f})"
          f" -> FWHM {r['T0']*1.665e15:.1f} fs")
    bw = 0.44 / (r['T0'] * 1.665)  # Gaussian TBP with FWHM
    print(f"  transform-limited optical bandwidth ~ {bw/1e12:.1f} THz")
