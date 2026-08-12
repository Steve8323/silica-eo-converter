#!/usr/bin/env python
"""Generate all report figures from out/*.json into fig/*.png."""
import glob
import json
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

os.makedirs("fig", exist_ok=True)
FW = dict(figsize=(7.2, 4.2), dpi=140)

def J(path):
    with open(path) as f:
        return json.load(f)

# ---------- 1. dispersion sweep ----------
if os.path.exists("out/dispersion.json"):
    d = J("out/dispersion.json")
    rows = d["rows"]
    w = [r["w"] for r in rows]
    b2 = [r["beta2"] * 1e27 for r in rows]
    T0 = [r["T0"] * 1.665e15 for r in rows]
    fig, ax = plt.subplots(**FW)
    ax.plot(w, b2, "-", color="#1f77b4")
    ax.axhline(0, color="k", lw=0.6)
    ax.axvline(d["best"]["w"], color="#d62728", ls="--", lw=1,
               label=f"zero-GVD w*={d['best']['w']:.3f} um")
    ax.axvline(0.744, color="#7f7f7f", ls=":", lw=1, label="single-mode limit 0.744 um")
    ax.set_xlabel("core width w (um)")
    ax.set_ylabel("beta2 (ps$^2$/km) @1550 nm")
    ax.set_title("Silica-slab waveguide GVD vs width (exact, Malitson Sellmeier)")
    ax2 = ax.twinx()
    ax2.plot(w, T0, "-", color="#2ca02c", alpha=0.7)
    ax2.set_ylabel("shortest 5 mm pulse FWHM (fs)", color="#2ca02c")
    ax2.set_ylim(0, 120)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig("fig/dispersion.png"); plt.close(fig)

# ---------- 2. FDTD dispersion validation ----------
if os.path.exists("out/meep_beta_check.json"):
    m = J("out/meep_beta_check.json")
    from dispersion import beta_derivs, slab_beta_TE0, C
    nus = np.array(m["nu"]); ng_m = np.array(m["ng"]); b2_m = np.array(m["beta2"])
    lam_um = C / nus * 1e6
    # analytic curves on same grid
    n1 = 1.444024 if m.get("const") else None
    b_an = np.array([slab_beta_TE0(l, 0.991, n1=n1) for l in lam_um])
    om = 2 * np.pi * nus
    p = np.polyfit(om - om.mean(), b_an, 7)
    ng_a = np.polyval(np.polyder(p, 1), om - om.mean()) * C
    b2_a = np.polyval(np.polyder(p, 2), om - om.mean())
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8), dpi=140)
    axs[0].plot(nus / 1e12, ng_m, label="MEEP FDTD")
    axs[0].plot(nus / 1e12, ng_a, "--", label="analytic slab")
    axs[0].set_xlabel("frequency (THz)"); axs[0].set_ylabel("group index")
    axs[0].legend(fontsize=8)
    axs[1].plot(nus / 1e12, b2_m * 1e27, label="MEEP FDTD")
    axs[1].plot(nus / 1e12, b2_a * 1e27, "--", label="analytic slab")
    axs[1].set_xlabel("frequency (THz)"); axs[1].set_ylabel("beta2 (ps$^2$/km)")
    axs[1].legend(fontsize=8)
    fig.suptitle("FDTD vs analytic slab dispersion (w=0.991 um, constant n; "
                 "offset = FDTD grid dispersion at res=25)")
    fig.tight_layout(); fig.savefig("fig/beta_check.png"); plt.close(fig)

# ---------- 3. 5 mm pulse propagation ----------
if os.path.exists("out/pulse_5mm.json"):
    from dispersion import slab_beta_TE0, C
    p5 = J("out/pulse_5mm.json")
    L = 5e-3; LAM0 = 1.55e-6; NU0 = C / LAM0
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 3.8), dpi=140)
    for ax, (wkey, style) in zip(axs, [("0.991", "#1f77b4"), ("0.744", "#ff7f0e")]):
        T0 = p5[wkey]["T0"]
        nus = NU0 + np.linspace(-60, 60, 4001) * 1e12
        lam = C / nus * 1e6
        b = np.array([slab_beta_TE0(l, float(wkey)) for l in lam])
        ok = ~np.isnan(b)
        nus2, b2v = nus[ok], b[ok]
        om = 2 * np.pi * nus2; om0 = 2 * np.pi * NU0
        Ein = np.exp(-T0 ** 2 * (om - om0) ** 2 / 2)
        pfit = np.polyfit(om - om0, b2v, 1)
        phi = (b2v - np.polyval(pfit, om - om0)) * L
        dom = om[1] - om[0]
        t = np.fft.fftshift(np.fft.fftfreq(len(om), d=dom / (2 * np.pi))) * 1e15
        Ii = np.abs(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Ein)))) ** 2
        Io = np.abs(np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(Ein * np.exp(1j * phi))))) ** 2
        ax.plot(t, Ii / Ii.max(), label="input")
        ax.plot(t, Io / Io.max(), label="after 5 mm")
        ax.set_xlim(-150, 150); ax.set_xlabel("time (fs)")
        ax.set_title(f"w={wkey} um, FWHM {p5[wkey]['fwhm']*1e15:.1f} fs, "
                     f"BW {p5[wkey]['bw_thz']:.1f} THz", fontsize=9)
        ax.legend(fontsize=8)
    axs[0].set_ylabel("normalized |E|^2")
    fig.suptitle("Shortest pulse surviving 5 mm (exact beta(omega) propagation)")
    fig.tight_layout(); fig.savefig("fig/pulse5mm.png"); plt.close(fig)

# ---------- 4. coupler design ----------
if os.path.exists("out/coupler_design.json"):
    rows = J("out/coupler_design.json")
    g = [r["gap"] for r in rows]; L50 = [r["L50"] for r in rows]
    rmin = [r["rmin"] for r in rows]; rmax = [r["rmax"] for r in rows]
    fig, ax = plt.subplots(**FW)
    ax.semilogy(g, L50, "o-", color="#1f77b4")
    ax.set_xlabel("coupler gap (um)"); ax.set_ylabel("50/50 length (um)", color="#1f77b4")
    ax2 = ax.twinx()
    ax2.fill_between(g, rmin, rmax, color="#2ca02c", alpha=0.25)
    ax2.axhline(0.5, color="#2ca02c", ls="--", lw=0.8)
    ax2.set_ylabel("splitting ratio over 26.8 THz band", color="#2ca02c")
    ax.set_title("Directional coupler: size & bandwidth vs gap (supermode, exact)")
    fig.tight_layout(); fig.savefig("fig/coupler.png"); plt.close(fig)

# ---------- 5. ring Q vs radius ----------
qfiles = [f for f in sorted(glob.glob("out/ring_q_R*.json")) if "_dn" not in f]
if qfiles:
    Rs, Qs = [], []
    for qf in qfiles:
        d = J(qf)
        if d["modes"]:
            m = max(d["modes"], key=lambda r: r["amp"])
            Rs.append(d["R"]); Qs.append(m["Q"])
    idx = np.argsort(Rs)
    fig, ax = plt.subplots(**FW)
    ax.semilogy(np.array(Rs)[idx], np.array(Qs)[idx], "o-")
    ax.set_xlabel("ring radius (um)"); ax.set_ylabel("bend-limited Q (FDTD Harminv)")
    ax.axhline(4e5, color="#d62728", ls="--", lw=1,
               label="Q needed @ 100 V/um Kerr shift")
    ax.set_title("Silica/air ring: radiation Q vs radius")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("fig/ring_q.png"); plt.close(fig)

# ---------- 6. ring resonance shift vs Delta-n (Harminv) ----------
sfiles = sorted(glob.glob("out/ring_q_R8_dn*.json")) + ["out/ring_q_R8.json"]
sfiles = [f for f in sfiles if os.path.exists(f)]
if len(sfiles) >= 2:
    pts = []
    for sf in sfiles:
        d = J(sf)
        dn = 0.0
        if "_dn" in sf:
            dn = float(sf.split("_dn")[1].replace(".json", ""))
        if d["modes"]:
            m = max(d["modes"], key=lambda r: r["amp"])
            pts.append((dn, m["lam_nm"], m["Q"]))
    pts.sort()
    dns = np.array([p[0] for p in pts]); lams = np.array([p[1] for p in pts])
    fig, ax = plt.subplots(**FW)
    ax.plot(dns * 1e3, lams - lams[0], "o-", label="FDTD (Harminv)")
    # analytic: dlam = lam * Gamma * dn / ng ; fit Gamma
    if len(dns) > 1 and dns[-1] > 0:
        G = np.polyfit(dns, lams - lams[0], 1)[0] / (lams[0] / 1.4948)
        dd = np.linspace(0, dns.max(), 50)
        ax.plot(dd * 1e3, lams[0] * G * dd / 1.4948, "--",
                label=f"dlam = lam*Gamma*dn/ng, Gamma={G:.2f}")
    ax.set_xlabel("electrode Delta-n (x1e-3)")
    ax.set_ylabel("resonance shift (nm)")
    ax.set_title("Master MRR: resonance shift vs electrode Delta-n (R=8 um, FDTD)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("fig/ring_shift.png"); plt.close(fig)

# ---------- 7. MZI transfer curve ----------
mzfiles = sorted(glob.glob("out/mzi_*_LC*.json"))
mz = [J(f) for f in mzfiles if f.endswith("_narrow.json")]
if mz:
    mz.sort(key=lambda d: d["dn"])
    dns = [d["dn"] for d in mz]
    lam = np.array(mz[0]["lam_nm"])
    i0 = np.argmin(np.abs(lam - 1550))
    s = [d["s"][i0] for d in mz]
    fig, ax = plt.subplots(**FW)
    ax.plot(dns, s, "o-", color="#1f77b4", label="FDTD")
    # analytic cos curve fit: s = -cos(2 pi Gamma dn Lp / lam)
    dd = np.linspace(0, max(dns), 200)
    from scipy.optimize import curve_fit
    def model(dn, G, s0):  # phase per dn
        return -np.cos(np.pi * dn / G + s0)
    try:
        popt, _ = curve_fit(model, dns, s, p0=[0.04, 0.0])
        ax.plot(dd, model(dd, *popt), "--", color="#d62728",
                label=f"-cos fit: dn_pi={popt[0]:.4f}")
    except Exception:
        pass
    ax.set_xlabel("electrode Delta-n in one arm")
    ax.set_ylabel("signed output s = (P_bar - P_x)/(P_bar + P_x)")
    ax.set_ylim(-1.05, 1.05)
    ax.set_title("MZI signed transfer curve @1550 (FDTD, couplers only)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("fig/mzi_transfer.png"); plt.close(fig)
    # broadband s(nu) if present
    bb = [J(f) for f in mzfiles if f.endswith("_broadband.json")]
    if bb:
        fig, ax = plt.subplots(**FW)
        for d in bb:
            lam = np.array(d["lam_nm"])
            ax.plot(lam, d["s"], label=f"dn={d['dn']:g}")
        ax.set_xlabel("wavelength (nm)"); ax.set_ylabel("signed output s")
        ax.set_title("MZI output vs wavelength (broadband FDTD)")
        ax.legend(fontsize=8)
        fig.tight_layout(); fig.savefig("fig/mzi_broadband.png"); plt.close(fig)

# ---------- 8. comb / FWM ----------
cfiles = sorted(glob.glob("out/mrr_comb_*.json"))
if cfiles:
    fig, ax = plt.subplots(**FW)
    for cf in cfiles:
        d = J(cf)
        fr = np.array(d["freq"]); psd = np.array(d["psd"])
        ax.semilogy(fr * 299.792458, psd / psd.max(), lw=0.8,
                    label=f"chi3={d['chi3']:g}")
        for fx, nm in [(d["F1"], "p1"), (d["F2"], "p2"),
                       (2 * d["F1"] - d["F2"], "2f1-f2"),
                       (2 * d["F2"] - d["F1"], "2f2-f1")]:
            ax.axvline(fx * 299.792458, color="k", lw=0.4, alpha=0.4)
    ax.set_xlabel("frequency (THz)"); ax.set_ylabel("PSD (norm.)")
    ax.set_title("Master MRR chi3 four-wave mixing: first comb sidebands (FDTD)")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig("fig/mrr_comb.png"); plt.close(fig)

print("figures:", sorted(os.listdir("fig")))
