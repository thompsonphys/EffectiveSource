"""Generate Taylor coefficients (in m1, complementary parameter) for
PK, QK, PE, QE where K(1-m1) = PK + QK*ln(1/m1), E(1-m1) = PE + QE*ln(1/m1),
QK = K(m1)/pi, QE = (K(m1)-E(m1))/pi.

Coefficients via Cauchy integral on |m1| = 0.5 (singularity at m1=1).
Also validate closed-form first/second m-derivative expressions vs FD.
"""
import mpmath as mp

mp.mp.dps = 60
N = 23          # highest kept order
R = mp.mpf("0.5")
NP = 512


def fset(m1):
    """Return PK, QK, PE, QE at complex m1."""
    Km = mp.ellipk(m1)
    Em = mp.ellipe(m1)
    Kc = mp.ellipk(1 - m1)
    Ec = mp.ellipe(1 - m1)
    QK = Km / mp.pi
    QE = (Km - Em) / mp.pi
    L = mp.log(1 / m1)
    return Kc - QK * L, QK, Ec - QE * L, QE


# Cauchy integral coefficients
vals = [fset(R * mp.exp(2j * mp.pi * k / NP)) for k in range(NP)]
coefs = {name: [] for name in ("pk", "qk", "pe", "qe")}
for idx, name in enumerate(("pk", "qk", "pe", "qe")):
    for n in range(N + 1):
        s = mp.mpc(0)
        for k in range(NP):
            z = R * mp.exp(2j * mp.pi * k / NP)
            s += vals[k][idx] / z**n
        c = s / NP
        assert abs(mp.im(c)) < mp.mpf("1e-40"), (name, n, c)
        coefs[name].append(mp.re(c))

# sanity: known values
assert abs(coefs["pk"][0] - 2 * mp.log(2)) < mp.mpf("1e-40")
assert abs(coefs["qk"][0] - mp.mpf(1) / 2) < mp.mpf("1e-40")
assert abs(coefs["pe"][0] - 1) < mp.mpf("1e-40")
assert abs(coefs["qe"][0]) < mp.mpf("1e-40")

# series accuracy at threshold m1=0.1 (value, f', f'')
for name in coefs:
    for m1 in (mp.mpf("0.02"), mp.mpf("0.1")):
        idx = ("pk", "qk", "pe", "qe").index(name)
        truth = fset(m1)[idx]
        val = sum(c * m1**n for n, c in enumerate(coefs[name]))
        d1t = mp.diff(lambda x: fset(x)[idx], m1)
        d1s = sum(n * c * m1 ** (n - 1) for n, c in enumerate(coefs[name]) if n)
        d2t = mp.diff(lambda x: fset(x)[idx], m1, 2)
        d2s = sum(n * (n - 1) * c * m1 ** (n - 2)
                  for n, c in enumerate(coefs[name]) if n > 1)
        print(f"{name} m1={float(m1):4} relerr val {float(abs(val/truth-1)):.1e}"
              f" d1 {float(abs(d1s/d1t-1)):.1e} d2 {float(abs(d2s/d2t-1)):.1e}")

# --- validate closed-form derivative expressions (the m1>=0.1 branch) --------
def closed(m):
    Km, Em = mp.ellipk(m), mp.ellipe(m)
    Kc, Ec = mp.ellipk(1 - m), mp.ellipe(1 - m)
    L = mp.log(1 / m)
    QK = Km / mp.pi
    QE = (Km - Em) / mp.pi
    PK = Kc - QK * L
    PE = Ec - QE * L
    Kmp = (Em - (1 - m) * Km) / (2 * m * (1 - m))
    Emp = (Em - Km) / (2 * m)
    Kcp = -(Ec - m * Kc) / (2 * m * (1 - m))
    Ecp = -(Ec - Kc) / (2 * (1 - m))
    QKp = Kmp / mp.pi
    QEp = (Kmp - Emp) / mp.pi
    PKp = Kcp - QKp * L + QK / m
    PEp = Ecp - QEp * L + QE / m
    # second derivatives
    D = 2 * m * (1 - m)
    Dp = 2 * (1 - 2 * m)
    Kmpp = ((Emp + Km - (1 - m) * Kmp) - Kmp * Dp) / D
    Empp = ((Emp - Kmp) - 2 * Emp) / (2 * m)
    Kcpp = ((-(Ecp - Kc - m * Kcp)) - Kcp * Dp) / D
    Ecpp = ((-(Ecp - Kcp)) - Ecp * (-2)) / (2 * (1 - m))
    QKpp = Kmpp / mp.pi
    QEpp = (Kmpp - Empp) / mp.pi
    PKpp = Kcpp - QKpp * L + 2 * QKp / m - QK / (m * m)
    PEpp = Ecpp - QEpp * L + 2 * QEp / m - QE / (m * m)
    return (PK, QK, PE, QE, PKp, QKp, PEp, QEp, PKpp, QKpp, PEpp, QEpp)


for m in (mp.mpf("0.1"), mp.mpf("0.3"), mp.mpf("0.7"), mp.mpf("0.95")):
    got = closed(m)
    worst = 0
    for idx in range(4):
        d1t = mp.diff(lambda x: fset(x)[idx], m)
        d2t = mp.diff(lambda x: fset(x)[idx], m, 2)
        worst = max(worst,
                    abs(got[4 + idx] / d1t - 1) if d1t else 0,
                    abs(got[8 + idx] / d2t - 1) if d2t else 0)
    print(f"closed-form derivs m1={float(m)}: worst rel {float(worst):.1e}")

# --- emit C arrays ------------------------------------------------------------
print("\n/* Taylor coefficients in m1 for PK, QK, PE, QE (|m1| < 0.1 branch) */")
print(f"#define EFFSOURCE_PQ_NSER {N + 1}")
for name in ("pk", "qk", "pe", "qe"):
    body = ",\n  ".join(mp.nstr(c, 26) for c in coefs[name])
    print(f"static const double effsource_{name}_ser[EFFSOURCE_PQ_NSER] = {{\n  {body}}};")
