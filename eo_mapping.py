#!/usr/bin/env python
"""Electrode physics: map simulated Delta-n values to real silica DC-Kerr.

Silica has no chi2 (centrosymmetric); with electrodes the only EO effect is
the DC Kerr effect:
    Delta n = 3 chi3 E^2 / (2 n)          chi3_1111 ~ 2.5e-22 m^2/V^2
Biased operation linearizes it:  E = Eb + e(t)  ->
    delta n(t) ~ (3 chi3 Eb / n) * e(t)   (effective Pockels)

Outputs the V -> Delta-n table, the ring Q needed to resolve each Delta-n,
the MZI L_pi, and the bandwidth budget of the converter.
"""
import json

CHI3 = 2.5e-22        # m^2/V^2, fused silica
N = 1.444024
NG = 1.4948
LAM = 1.55e-6
NU = 299792458.0 / LAM
GAMMA = 0.8           # field overlap of optical mode with electrode E region
E_BD = 1.0e9          # V/m, thin-film silica breakdown (~1 GV/m)

def dn_kerr(E):
    return 3 * CHI3 * E ** 2 / (2 * N)

rows = []
print(f"{'E [V/um]':>9} {'Dn_Kerr':>10} {'ring Q_min':>11} {'MZI L_pi':>10}")
for E_Vum in (10, 50, 100, 200, 300, 500):
    E = E_Vum * 1e6
    dn = dn_kerr(E)
    # ring: shift >= half linewidth  ->  Q >= ng lam? shift dnu/nu = -Gamma dn/ng
    #   dnu = nu Gamma dn / ng ; linewidth = nu/Q ; require shift >= linewidth/2
    Qmin = NG / (2 * GAMMA * dn)
    Lpi = LAM / (2 * GAMMA * dn)
    rows.append(dict(E_Vum=E_Vum, dn=dn, Qmin=Qmin, Lpi_m=Lpi))
    print(f"{E_Vum:9g} {dn:10.3e} {Qmin:11.3e} {Lpi*100:8.2f} cm")

print(f"\nbreakdown limit ~{E_BD/1e6:.0f} V/um -> max static Dn = {dn_kerr(E_BD):.2e}")
print("silica wedge/toroid resonators reach Q = 1e7-1e9 (Kerr-comb literature),")
print("so E ~ 50-100 V/um is resolvable by the master MRR.")

# photon-lifetime-limited electrical bandwidth vs Q
print(f"\n{'Q':>8} {'f_3dB = nu/2Q':>14}")
bw = {}
for Q in (1e5, 3e5, 1e6, 1e7, 1e8):
    f3 = NU / (2 * Q)
    bw[f"{Q:g}"] = f3
    print(f"{Q:8g} {f3/1e9:11.3f} GHz")

with open("out/eo_mapping.json", "w") as f:
    json.dump(dict(rows=rows, chi3=CHI3, gamma=GAMMA, bw=bw), f, indent=1)
