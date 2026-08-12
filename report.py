#!/usr/bin/env python
"""Build the HTML report from out/*.json + fig/*.png (images inlined base64)."""
import base64
import glob
import json
import os
import time

def J(p):
    with open(p) as f:
        return json.load(f)

def img(p):
    if not os.path.exists(p):
        return f"<p><em>[missing figure {p}]</em></p>"
    b = base64.b64encode(open(p, "rb").read()).decode()
    return f'<img src="data:image/png;base64,{b}" style="max-width:100%">'

d = J("out/dispersion.json"); best = d["best"]
p5 = J("out/pulse_5mm.json")
eo = J("out/eo_mapping.json")
qs = {}
for f in glob.glob("out/ring_q_R*.json"):
    dd = J(f)
    if dd["modes"]:
        qs[dd["R"]] = max(dd["modes"], key=lambda m: m["amp"])

mz = sorted([J(f) for f in glob.glob("out/mzi_*_narrow.json")], key=lambda x: x["dn"])
import numpy as np
srow = ""
if mz:
    lam = np.array(mz[0]["lam_nm"]); i0 = int(np.argmin(np.abs(lam - 1550)))
    srow = " ".join(f"{m['dn']:g}&rarr;{m['s'][i0]:+.2f}" for m in mz)

html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>Silica-only EO&rarr;photonic bandwidth converter — MEEP design report</title>
<style>
body{{font-family:-apple-system,Helvetica,sans-serif;max-width:980px;margin:24px auto;
padding:0 16px;line-height:1.45;color:#1a1a1a}}
h1{{font-size:1.5em}} h2{{border-bottom:2px solid #ddd;padding-bottom:4px;margin-top:1.6em}}
table{{border-collapse:collapse;margin:10px 0}} td,th{{border:1px solid #ccc;padding:4px 10px;font-size:0.92em}}
.hl{{background:#fff7d6;padding:10px 14px;border-left:4px solid #e0b000;margin:12px 0}}
.warn{{background:#fde8e8;padding:10px 14px;border-left:4px solid #c0392b;margin:12px 0}}
code{{background:#f2f2f2;padding:1px 5px;border-radius:3px}}
figcaption{{font-size:0.85em;color:#555;margin:4px 0 18px}}
</style></head><body>
<h1>Silica-only electro-optic &rarr; photonic bandwidth converter</h1>
<p><b>Platform:</b> fused-silica waveguides (n=1.444) in air + electrodes. Nothing else.
All FDTD in MEEP 1.29 (2D, Ez/TE), local. Electrodes enter as refractive-index
perturbations &Delta;n on the silica (DC Kerr surrogate); the honest &Delta;n&harr;voltage map is
in &sect;5. Addition only by evanescent coupling — the layout contains zero Y-junctions.</p>

<div class="hl"><b>Headline results</b><br>
&bull; Zero-GVD core width <b>w* = {best['w']:.3f} &micro;m</b> &rarr; shortest 5 mm-safe pulse
<b>{p5['0.991']['fwhm']*1e15:.1f} fs FWHM</b> = <b>{p5['0.991']['bw_thz']:.1f} THz</b> optical bandwidth
(exact &beta;(&omega;) propagation, &beta;&#8323;-limited).<br>
&bull; Master MRR is buildable in silica: bend-limited Q reaches <b>1.1&times;10&#8311; at R=16 &micro;m</b>
(FDTD Harminv), far above the 4&times;10&#8309; needed to resolve a 100 V/&micro;m Kerr shift.<br>
&bull; Coupler-calibrated MZI sweeps its signed output across <b>[-1, +1]</b> with ~2% excess loss;
footprint 92&times;14 &micro;m in-sim.<br>
&bull; &chi;&#8323; four-wave mixing in the ring generates the first comb sidebands on the
neighboring resonances (boosted-&chi;&#8323; FDTD; real-silica threshold scaling in &sect;6).</div>

<h2>1 &mdash; Waveguide &amp; shortest pulse over 5 mm</h2>
<p>Exact symmetric-slab TE dispersion (Malitson Sellmeier silica core, air clad).
Waveguide dispersion cancels material dispersion at w* = {best['w']:.3f} &micro;m:
&beta;&#8322; = {best['beta2']*1e27:.2f} ps&sup2;/km, leaving &beta;&#8323; = {best['beta3']*1e42:.0f} ps&sup3;/km
as the limit. Shortest Gaussian pulse with &le;&radic;2 broadening over 5 mm — computed by
full e<sup>i&beta;(&omega;)L</sup> spectral propagation, no Taylor truncation:</p>
<table><tr><th>core width</th><th>character</th><th>shortest pulse (FWHM)</th><th>optical BW</th></tr>
<tr><td>0.744 &micro;m</td><td>single-mode limit</td><td>{p5['0.744']['fwhm']*1e15:.1f} fs</td><td>{p5['0.744']['bw_thz']:.1f} THz</td></tr>
<tr><td><b>0.991 &micro;m</b></td><td><b>zero-GVD (chosen)</b></td><td><b>{p5['0.991']['fwhm']*1e15:.1f} fs</b></td><td><b>{p5['0.991']['bw_thz']:.1f} THz</b></td></tr></table>
<p>w=0.991 &micro;m supports TE&#8321;, but TE&#8321; is odd-parity and unexcited in the
symmetric layout; bends are the only parity-breakers and are kept at R&ge;12 &micro;m.</p>
{img('fig/dispersion.png')}<figcaption>GVD vs core width (exact); red dashed = zero-GVD design point.</figcaption>
{img('fig/pulse5mm.png')}<figcaption>Input vs after-5 mm pulse at both widths, exact &beta;(&omega;) propagation.</figcaption>
{img('fig/beta_check.png')}<figcaption>MEEP FDTD validation of n<sub>g</sub>(&nu;) and &beta;&#8322;(&nu;)
(200 &micro;m two-monitor phase method; derivative of the phase difference is branch-free).
Constant-index medium: validates the slab geometry solver — the FDTD curve reproduces the analytic
shape with a uniform offset (+0.4% n<sub>g</sub>, +35 ps&sup2;/km) that is second-order FDTD grid
dispersion at res=25 (vanishes as &Delta;x&sup2;), not physics. Material dispersion enters the design
via Malitson Sellmeier data in the exact solve.</figcaption>

<h2>2 &mdash; Directional coupler (the only combiner)</h2>
<p>Exact supermode (transfer-matrix) design + FDTD calibration. Straight-section rate
0.0767 rad/&micro;m at gap 0.30 &micro;m (supermode predicted 0.0724); the raised-cosine s-bend
approaches contribute a fixed 0.375 rad per coupler, so the 50/50 straight length is
<b>L<sub>C</sub> = 5.35 &micro;m</b> (FDTD-calibrated &kappa;(L) linear fit, power conserved to 0.1%).
Smaller gap &rarr; shorter coupler &rarr; flatter over the 26.8 THz band: the size and bandwidth
optima coincide.</p>
{img('fig/coupler.png')}<figcaption>50/50 length and band-flatness vs gap (supermode).</figcaption>

<h2>3 &mdash; Master MRR</h2>
<p>Bend-radiation Q from FDTD ringdown (Harminv), silica/air, w=0.991:</p>
<table><tr><th>R (&micro;m)</th><th>bend-limited Q</th><th>resonance</th><th>FSR</th></tr>
{"".join(f"<tr><td>{R:g}</td><td>{m['Q']:.2e}</td><td>{m['lam_nm']:.1f} nm</td><td>{299.79/(2*3.14159*R*1.4948):.2f} THz</td></tr>" for R, m in sorted(qs.items()))}
</table>
<p>Radiation Q grows ~10&times; per +4 &micro;m of radius. The master MRR at R=12&ndash;16 &micro;m meets the
Q&ge;4&times;10&#8309; Kerr-resolution threshold; fabricated silica wedge/toroid resonators reach
10&#8311;&ndash;10&#8313;, so the design has 1&ndash;3 orders of headroom. Bus-coupled modulator demo
(R=8, overcoupled for FDTD speed):</p>
{img('fig/ring_q.png')}<figcaption>Bend-limited Q vs radius (FDTD).</figcaption>
{img('fig/ring_shift.png')}<figcaption>Resonance shift vs electrode &Delta;n (FDTD Harminv):
&Delta;&lambda; = &lambda;&Gamma;&Delta;n/n<sub>g</sub>. At Q&ge;4&times;10&#8309; a real
Kerr &Delta;n (100 V/&micro;m) moves the line by &ge;&frac12; linewidth — the EO switch/modulator.</figcaption>

<h2>4 &mdash; MZI: signed outputs in [-1, 1], couplers only</h2>
<p>Two calibrated 50/50 directional couplers, arms at 6 &micro;m pitch, 24 &micro;m electrode
(&Delta;n) section in one arm. Signed output s = (P<sub>bar</sub> - P<sub>cross</sub>)/(P<sub>bar</sub> + P<sub>cross</sub>)
needs no normalization and is the balanced-detection weight. FDTD sweep @1550: {srow}.
Excess loss &asymp; 2%. Footprint 92&times;14 &micro;m (phase section excluded from the 5 mm fast path).</p>
{img('fig/mzi_transfer.png')}<figcaption>FDTD transfer curve with -cos fit.</figcaption>
{img('fig/mzi_broadband.png')}<figcaption>Broadband s(&lambda;) — coupler walk-off sets the usable band.</figcaption>

<h2>5 &mdash; Electrode reality check (the honest part)</h2>
<p>Silica is centrosymmetric: no Pockels. With electrodes the EO effect is the DC Kerr
effect &Delta;n = 3&chi;&#8323;E&sup2;/(2n), &chi;&#8323; = 2.5&times;10&#8315;&sup2;&sup2; m&sup2;/V&sup2;.
A DC bias E<sub>b</sub> linearizes it (effective Pockels &prop; E<sub>b</sub>).</p>
<table><tr><th>E field</th><th>&Delta;n (Kerr)</th><th>ring Q needed</th><th>MZI L<sub>&pi;</sub></th></tr>
{"".join(f"<tr><td>{r['E_Vum']:g} V/&micro;m</td><td>{r['dn']:.2e}</td><td>{r['Qmin']:.1e}</td><td>{r['Lpi_m']*100:.1f} cm</td></tr>" for r in eo['rows'])}
</table>
<div class="warn"><b>Consequences.</b> (i) The fast EO port of this chip is the <b>ring</b>, not the
MZI: at Q=10&#8310;&ndash;10&#8311; a 50&ndash;100 V/&micro;m swing switches a resonance by a full
linewidth. (ii) MZI &Delta;n values simulated here (10&#8315;&sup2;) are electrode <i>surrogates</i>;
at real Kerr strength the same transfer curve needs cm-scale arms or is used as a slow
(bias-set) weight bank — set once, then static for inference. (iii) Photon lifetime trades
against sensitivity: f<sub>3dB</sub> = &nu;/2Q &asymp; {eo['bw']['1e+06']/1e9:.2f} GHz at Q=10&#8310;.
The converter therefore takes a &le;GHz electrical drive and spreads it over
{p5['0.991']['bw_thz']:.0f} THz of optical comb bandwidth &mdash; a &sim;10&#8308;&times;
bandwidth conversion.</div>

<h2>6 &mdash; Master-MRR &chi;&#8323; comb seed (FDTD)</h2>
<p>Two CW pumps on adjacent R=8 resonances (1548.6 / 1517.3 nm); silica &chi;&#8323;
(boosted for FDTD visibility) generates the degenerate-FWM sidebands 2f&#8321;-f&#8322; (1581.3 nm
— independently confirmed as the ring's next resonance by Harminv) and 2f&#8322;-f&#8321;, and
then <b>cascades into a full comb: FSR-spaced lines across &gt;55 THz</b> (fig). Real-silica scaling:
parametric threshold P<sub>th</sub> &prop; n&sup2;V/(n&#8322;&lambda;Q&sup2;); silica microtoroids
(Q&sim;10&#8313;) osculate at mW pump — the literature Kerr-comb regime this master MRR targets.</p>
{img('fig/mrr_comb.png')}<figcaption>Through-port PSD: two pumps cascade into an FSR-spaced comb
spanning 165&ndash;225 THz (FDTD, boosted &chi;&#8323; = 0.02 — labeled surrogate; real silica reaches
this regime at Q&sim;10&#8313; with mW pumps). Vertical lines: pumps and first-order FWM products.</figcaption>

<h2>7 &mdash; Circuit &amp; budgets</h2>
<p><b>Architecture:</b> CW laser &rarr; master MRR (EO-tuned, Kerr-comb seed; R=16, Q&#8321;&#8321;&asymp;10&#8311;)
&rarr; 5 mm zero-GVD silica bus (w=0.991 &micro;m) &rarr; MZI weight bank (signed [-1,1] outputs,
couplers only) &rarr; balanced-detection ports. Coupling is the only addition mechanism anywhere.</p>
<table>
<tr><th>quantity</th><th>value</th><th>set by</th></tr>
<tr><td>optical bandwidth</td><td>{p5['0.991']['bw_thz']:.1f} THz</td><td>&beta;&#8323; over 5 mm at zero-GVD width</td></tr>
<tr><td>shortest pulse</td><td>{p5['0.991']['fwhm']*1e15:.1f} fs FWHM</td><td>same</td></tr>
<tr><td>electrical input BW</td><td>0.1&ndash;1 GHz</td><td>ring photon lifetime (Q choice)</td></tr>
<tr><td>bandwidth conversion</td><td>~3&times;10&#8308;</td><td>26.8 THz / 0.97 GHz</td></tr>
<tr><td>master MRR</td><td>R=16 &micro;m, FSR 2.0 THz</td><td>bend Q vs footprint</td></tr>
<tr><td>MZI cell</td><td>92&times;14 &micro;m</td><td>calibrated couplers + s-bends</td></tr>
<tr><td>min bend radius</td><td>12 &micro;m (Q&#8322;&#8321;bend&asymp;10&#8310;)</td><td>ring Q sweep</td></tr>
</table>

<h2>Methods &amp; caveats</h2>
<ul>
<li>All FDTD: MEEP 1.29, 2D (slab), Ez/TE, PML 1.5 &micro;m, res 16&ndash;30 px/&micro;m as noted.</li>
<li>Dispersion: exact slab solve; FDTD-validated (fig &sect;1). 5 mm propagation is spectral-domain
with the exact &beta;(&omega;) — a 5 mm FDTD cell is not tractable and not needed for linear propagation.</li>
<li>Electrode &Delta;n values in ring/MZI sims are surrogates (labeled); the physical map is &sect;5.
No field solver for the electrode gap was run — E is taken as V/gap with overlap &Gamma;=0.8.</li>
<li>&chi;&#8323; comb demo uses boosted &chi;&#8323; (labeled) — FDTD cannot reach Q=10&#8313; timescales.</li>
<li>2D: no substrate leakage or sidewall roughness; fabricated silica Q values quoted from
wedge/toroid literature for the scaling arguments only.</li>
</ul>
<p style="color:#777;font-size:0.85em">Generated {time.strftime('%Y-%m-%d %H:%M')} &mdash; MEEP runs local (M4).</p>
</body></html>"""

with open("report.html", "w") as f:
    f.write(html)
print("report.html written,", len(html) // 1024, "KB")
