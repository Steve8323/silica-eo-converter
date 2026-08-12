# Silica-only EO → photonic bandwidth converter

Fused-silica waveguides in air + electrodes — nothing else. All FDTD in MEEP
1.29 (2D Ez/TE, local `mp` conda env). Electrodes enter the sims as Δn
surrogates on the silica; the honest Δn↔voltage map (silica DC Kerr,
χ³ = 2.5e-22 m²/V²) is `eo_mapping.py` / report §5.

## Headline results

| quantity | value | from |
|---|---|---|
| zero-GVD core width | **w\* = 0.991 µm** | exact slab solve (`dispersion.py`) |
| shortest 5 mm-safe pulse | **16.4 fs FWHM = 26.8 THz** | exact e^{iβ(ω)L} propagation (`pulse_5mm.py`) |
| ring bend Q | 2.6e3 (R=5) → **1.1e7 (R=16)** | FDTD Harminv (`ring_q.py`) |
| ring EO tuning | Δλ = λΓΔn/n_g, **Γ = 0.97**, linear | FDTD (`ring_q.py R dn`) |
| 50/50 coupler | gap 0.3 µm, **L_C = 5.35 µm** (s-bends included) | FDTD calibration (`mzi.py 0 L coupler`) |
| MZI signed output | **s ∈ [−0.98, +0.97]**, Δn_π = 0.0332, ~2% loss | FDTD sweep (`mzi.py`) |
| χ³ comb | 2 pumps → FSR-spaced comb over >55 THz | nonlinear FDTD, boosted χ³ (`mrr_comb.py`) |
| bandwidth conversion | ≤1 GHz electrical → 26.8 THz optical (~3e4×) | photon lifetime vs pulse BW |

Architecture: CW laser → master MRR (R=16, EO-tuned, Kerr-comb seed) → 5 mm
zero-GVD bus → MZI weight bank (signed [−1,1], directional couplers only — no
waveguide merges anywhere) → balanced detection.

## Reproduce

```bash
PY=/opt/anaconda3/envs/mp/bin/python
$PY dispersion.py            # exact slab GVD sweep -> out/dispersion.json
$PY pulse_5mm.py             # shortest pulse over 5 mm
$PY meep_beta_check.py const # FDTD validation of the slab solver
$PY coupler_design.py        # supermode coupler table
$PY ring_q.py 8 [dn]         # ring Q / resonance shift
$PY mzi.py DN 5.35 [coupler|broadband]
$PY mrr_comb.py 0.64574 0.65909 0.02 1500
$PY figs.py && $PY report.py # figures + report.html
```

## Caveats (see report §Methods)

- 2D slab physics throughout; no substrate leakage / roughness.
- 5 mm propagation is exact spectral-domain with FDTD-validated β(ω) — a 5 mm
  FDTD cell is intractable and unnecessary for linear propagation.
- Ring/MZI Δn values are electrode surrogates (real Kerr Δn ~ 2.6e-6 at
  100 V/µm → the fast EO port is the high-Q ring; MZI weights are slow/bias-set
  or cm-scale).
- χ³ in `mrr_comb.py` is boosted for FDTD timescales (labeled); real-silica
  comb operation needs the literature Q~1e8-1e9 wedge/toroid regime.
