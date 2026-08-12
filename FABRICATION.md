# Fabrication flow — N=50 silica-only EO→photonic converter

Platform: **suspended thermal-oxide membrane on Si** — the waveguide core is
undoped SiO₂ with air on all sides (this is what makes it "silica + electrodes,
nothing else"; there is no lower-index solid to clad silica with). Three mask
levels. Minimum feature: **0.12 µm coupler gap** → e-beam or DUV (193 nm)
lithography for mask 1; masks 2–3 are ≥1 µm and can be stepper/contact.

## Mask set (artwork: `out/masks.json`, previews `fig/mask_*.png`, interactive `viewer.html`)

| # | Name | Field | Purpose | Min feature |
|---|---|---|---|---|
| 1 | SEO-WG-01 | dark | waveguide etch | 0.12 µm gap |
| 2 | SEO-MET-02 | clear | electrode liftoff | 1.0 µm spacing |
| 3 | SEO-REL-03 | dark | XeF₂ release windows | 6 µm |

All three carry corner cross alignment marks; overlay budget 0.25 µm
(electrode-to-guide spacing tolerates ±0.5 µm; coupler gaps are defined
entirely within mask 1, so overlay does NOT touch the critical dimensions).

## Process steps

1. **Substrate**: Si ⟨100⟩ wafer. Grow **1.0 µm thermal oxide** (wet, 1100 °C).
   This oxide IS the device layer — its thickness sets the out-of-plane mode.
2. **Mask 1 — waveguides (SEO-WG-01)**: spin e-beam resist (or DUV stack);
   expose buses, master MRR (R=16, ring–bus gap 0.15), all 50 ladder cells with
   the per-cell gaps/lengths from `out/pattern50.json`, facet inverse tapers,
   dumps. **CD-critical**: tap gaps 0.12–0.5 µm, coupler lengths ±50 nm.
3. **Oxide etch**: RIE full etch through the 1.0 µm oxide (C₄F₈/O₂ or
   CHF₃/Ar), vertical sidewalls. Strip resist.
4. **Sidewall smoothing anneal**: 1000 °C N₂, 2 h (reflow-smooths etch
   roughness — silica's low index makes it tolerant, but ring Q benefits).
5. **Mask 2 — electrodes (SEO-MET-02)**: image-reversal resist, liftoff of
   **Ti 5 nm / Au 150 nm**: MRR inner+outer ring electrodes (3 µm wide, 1 µm
   from the guide edge), RF + DC pads, 50 weight-electrode pairs flanking each
   arm's vertical run, signal/ground rails and fan-out. Metal never crosses a
   suspended span (rails run on unreleased mesa).
6. **Mask 3 — release (SEO-REL-03)**: pattern release windows (resist as
   mask), then **XeF₂ isotropic Si etch** undercuts ~8 µm per side: buses,
   couplers, ring and arms end up as an air-clad suspended silica membrane
   held by the unetched field + window spacing (24 µm pitch → continuous
   tether ribs). Resist strip in O₂ plasma (no liquids after release).
7. **Facets**: deep-etched facet trench is included in mask 1 at both left
   facets (input bus A taper, output bus B taper); cleave through the trench.
   Inverse tapers (0.3 µm tip) mode-match to lensed fiber.
8. **Anneal + CPO**: final 400 °C bake to drive off water (silica OH loss).

## Electrical hookup

- **RF pad** (master MRR outer ring) + **DC bias pad** (inner ring): bias
  E_b ≈ 50–100 V/µm across the 3 µm gap linearizes the DC Kerr effect; the
  GHz signal rides on top. Modulation bandwidth = photon lifetime limit
  (ν/2Q ≈ 0.1–1 GHz at the working Q).
- **Weight pads** (rail + per-arm fingers): quasi-static, set once per weight
  configuration. Also used for the flatness trim (below).

## Post-fab calibration (mirrors what the FDTD needed)

1. Broadband source in, camera on the dumps: measure per-cell tap ratio;
   compare to the 1/(51−k) targets.
2. Two trim passes on the weight electrodes (the N=6 FDTD went ±3.3 dB →
   ±1.5 dB in exactly two passes; electrical trim replaces geometric trim).
3. Lock the master MRR to the laser with the DC pad (thermal drift of silica
   is small: dn/dT ≈ 1e-5/K, but Q=10⁶ lines are 0.2 pm — budget a slow lock).

## Known risks

- 0.12 µm gaps at 1 µm depth: aspect ratio 8 in the gap — RIE lag; bias the
  mask CD −20 nm and verify on SEM cross-sections.
- Released membrane stress: thermal oxide is compressive (~300 MPa) —
  the 24 µm window pitch keeps spans short; buckling shows up as slot-time
  jitter first (measure slot times vs design as the stress metric).
- Au near the mode: 1 µm edge spacing costs <0.01 dB/cm at 1550 (evanescent
  tail is ~200 nm); do not shrink below 0.8 µm.
