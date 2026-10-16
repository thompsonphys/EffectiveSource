"""Check the analytic m-modes written by mmode-check --dump FILE against mpmath: mmode_exact_check.py FILE."""

import sys

import mpmath as mp

MMAX = 20
MAX_ERROR = 1e-15
NEIGHBOUR_FLOOR = 1e-2
FOURIER_C1 = 0.03


def kernel(m, nu, A, B):
    """K_m^(nu) = int_0^2pi cos(m psi) (A - B cos psi)^-nu dpsi from the hypergeometric function."""
    return (2 * mp.pi * A ** (-nu) * (B / (2 * A)) ** m * mp.rf(nu, m) / mp.factorial(m)
            * mp.hyp2f1((nu + m) / 2, (nu + m + 1) / 2, m + 1, (B / A) ** 2))


def to_s2_basis(coefficients, alpha, beta):
    """Coefficients of sum_j c_j u^j in powers of s^2 = alpha + beta u."""
    out = [mp.mpf(0)] * len(coefficients)
    for j, c in enumerate(coefficients):
        for k in range(j + 1):
            out[k] += c * mp.binomial(j, k) * (-alpha) ** (j - k) / beta ** j
    return out


def u_power_fourier(j):
    """Complex Fourier coefficients of u^j = ((1 - cos psi) / 2)^j."""
    coefficients = {0: mp.mpf(1)}
    for _ in range(j):
        product = {}
        for l, c in coefficients.items():
            for dl, w in ((0, mp.mpf(1) / 2), (1, -mp.mpf(1) / 4), (-1, -mp.mpf(1) / 4)):
                product[l + dl] = product.get(l + dl, 0) + c * w
        coefficients = product
    return coefficients


def exact_modes(alpha, beta, phi0, P, Q):
    """S_m for m = 0..MMAX in the better conditioned basis, at the working precision."""
    A, B = alpha + beta / 2, beta / 2
    if alpha / beta < FOURIER_C1:
        nus = [mp.mpf(11) / 2 - k for k in range(9)]
        K = {(m, k): kernel(m, nus[k], A, B) for m in range(MMAX + 1) for k in range(9)}
        e = to_s2_basis(P, alpha, beta)
        G = [mp.mpf(0)] * 9
        for j, q in enumerate(Q):
            G[j + 1] += 2 * q
            G[j + 2] -= 2 * q
        g = to_s2_basis(G, alpha, beta)
        modes = []
        for m in range(MMAX + 1):
            even = mp.fsum(e[k] * K[(m, k)] for k in range(len(e)))
            odd = mp.fsum(g[k] * mp.fsum(K[(abs(m - 1 - 2 * l), k)] for l in range(m)) for k in range(len(g)))
            modes.append(mp.expj(-m * phi0) * (even - 1j * odd))
        return modes
    K = [kernel(m, mp.mpf(11) / 2, A, B) for m in range(MMAX + 9)]
    numerator = {}
    for j, p in enumerate(P):
        for l, c in u_power_fourier(j).items():
            numerator[l] = numerator.get(l, 0) + p * c
    for j, q in enumerate(Q):
        for l, c in u_power_fourier(j).items():
            for dl, w in ((1, 1 / mp.mpc(0, 4)), (-1, -1 / mp.mpc(0, 4))):
                numerator[l + dl] = numerator.get(l + dl, 0) + q * c * w
    return [mp.expj(-m * phi0) * mp.fsum(c * K[abs(m - l)] for l, c in numerator.items()) for m in range(MMAX + 1)]


def local_scale(modes, m):
    """|S_m| plus NEIGHBOUR_FLOOR times its neighbours, so that an accidental zero of one mode is not graded relative to itself."""
    neighbours = [abs(modes[k]) for k in (m - 1, m + 1) if 0 <= k < len(modes)]
    return abs(modes[m]) + NEIGHBOUR_FLOOR * max(neighbours)


def exact_double(text):
    """The exact value of a C99 hexadecimal double."""
    return mp.mpf(float.fromhex(text))


def parse(path):
    """Cases (tag, alpha, beta, phi0, P, Q, S) of a dump file."""
    lines = open(path).read().splitlines()
    cases = []
    for i in range(0, len(lines), 4):
        head, p, q, s = (line.split() for line in lines[i:i + 4])
        pairs = lambda f: [exact_double(f[1 + 2 * j]) + exact_double(f[2 + 2 * j]) for j in range((len(f) - 1) // 2)]
        modes = [mp.mpc(exact_double(s[1 + 2 * m]), exact_double(s[2 + 2 * m])) for m in range(MMAX + 1)]
        cases.append((head[0], exact_double(head[1]) + exact_double(head[2]), exact_double(head[3]),
                      exact_double(head[4]), pairs(p), pairs(q), modes))
    return cases


def main():
    """Print the worst relative error of every mode per decade of C1 and exit with status 1 above MAX_ERROR."""
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    mp.mp.dps = 45
    worst = {}
    for tag, alpha, beta, phi0, P, Q, modes in parse(sys.argv[1]):
        exact = exact_modes(alpha, beta, phi0, P, Q)
        error = max(float(abs(modes[m] - exact[m]) / local_scale(exact, m)) for m in range(MMAX + 1))
        decade = int(mp.floor(mp.log10(alpha / beta)))
        worst[(tag, decade)] = max(worst.get((tag, decade), 0.0), error)
    print(f"analytic m-modes m <= {MMAX} vs mpmath, worst error of a mode relative to |S_m| + {NEIGHBOUR_FLOOR} max(|S_m-1|, |S_m+1|)")
    for tag, decade in sorted(worst):
        print(f"  {tag:3s} C1 in [1e{decade}, 1e{decade + 1}): {worst[(tag, decade)]:.1e}")
    overall = max(worst.values())
    print(f"check exact_modes {'PASS' if overall <= MAX_ERROR else 'FAIL'} {overall:.2e}")
    sys.exit(0 if overall <= MAX_ERROR else 1)


if __name__ == "__main__":
    main()
