#!/usr/bin/env python
"""Exact 1D multilayer TE mode solver (transfer matrix) for 2D-slab stacks.

Stack is a list of (n, thickness_um) with semi-infinite outer claddings given
separately.  Finds guided TE (Ez) modes: field ~ exp(-gamma |y|) outside,
transfer matrix through the layers, root-find the dispersion determinant.
Used for the directional-coupler supermodes (air|core|gap|core|air).
"""
import numpy as np

def _det(neff, k0, n_out, layers):
    """Determinant whose zeros are guided modes (decaying both sides)."""
    g = k0 * np.sqrt(neff ** 2 - n_out ** 2 + 0j)
    # start with decaying solution on the left: E = exp(+g y) -> [E, E'] = [1, g]
    v = np.array([1.0 + 0j, g])
    for n, d in layers:
        k = k0 * np.sqrt(n ** 2 - neff ** 2 + 0j)
        c, s = np.cos(k * d), np.sin(k * d)
        M = np.array([[c, s / k], [-k * s, c]]) if abs(k) > 1e-12 else \
            np.array([[1, d], [0, 1]])
        v = M @ v
    # right side must be purely decaying: E' = -g E
    return np.real(v[1] + g * v[0])

def modes(k0, n_out, layers, n_max, n_pts=4000):
    """All guided TE mode neffs in (n_out, n_max), descending."""
    lo, hi = n_out + 1e-6, n_max - 1e-6
    xs = np.linspace(lo, hi, n_pts)
    fs = np.array([_det(x, k0, n_out, layers) for x in xs])
    roots = []
    from scipy.optimize import brentq
    for i in range(len(xs) - 1):
        if np.isfinite(fs[i]) and np.isfinite(fs[i + 1]) and fs[i] * fs[i + 1] < 0:
            try:
                r = brentq(lambda x: _det(x, k0, n_out, layers), xs[i], xs[i + 1],
                           xtol=1e-13)
                roots.append(r)
            except ValueError:
                pass
    return sorted(roots, reverse=True)

def coupler_supermodes(lam_um, w, gap, n_core):
    """(even, odd) supermode neff of two identical slabs separated by gap."""
    k0 = 2 * np.pi / lam_um
    layers = [(n_core, w), (1.0, gap), (n_core, w)]
    ms = modes(k0, 1.0, layers, n_core)
    if len(ms) < 2:
        return None
    return ms[0], ms[1]

def coupling_length(lam_um, w, gap, n_core):
    """Full-transfer length L_pi = lam / (2 dn); 50/50 at L_pi/2."""
    sm = coupler_supermodes(lam_um, w, gap, n_core)
    if sm is None:
        return None
    ne, no = sm
    return lam_um / (2 * (ne - no))
