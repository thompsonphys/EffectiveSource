"""Replace the double-precision calc_m_split machinery with a long-double
core.  The quotient-rule assembly cancels ~1/C1 near the particle; long
double buys 2200x on that noise floor (srcA ~1e-14/C1).

Strips the previously spliced regions, then re-splices:
  - ld series coefficients (26 digits, L literals)
  - effsource_ellPQ_derivs_from_C1l (ld) + double wrapper (for tests)
  - ld box helper, ld prelude/den (type-transformed verbatim clones),
    ld num loops and channel assembly
"""
import sys, os

CODEGEN = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, CODEGEN)

PATH = "/home/jonathan/projects/sf/effectivesource/kerr-equatorial-ctx.c"
lines = open(PATH).read().split("\n")


def find_line(pat, start=0):
    for i in range(start, len(lines)):
        if pat in lines[i]:
            return i
    raise SystemExit(f"not found: {pat}")


def close_brace(i):
    while lines[i] != "}":
        i += 1
    return i


# --- strip old regions (delete higher-index region first) --------------------
bA = find_line("/* Taylor coefficients in m1")
# delete through the DOUBLE wrapper's close (its "(double C1" signature comes
# after the long-double "...C1l" primitive); both live in the spliced block.
eA = close_brace(find_line("void effsource_ellPQ_derivs_from_C1(double C1", bA))
bB = find_line("/* Scalar wave operator on one channel")
eB = close_brace(find_line("void effsource_equatorial_ctx_calc_m_split(", bB))
assert bB > eA
del lines[bB:eB + 1]
del lines[bA:eA + 1]
src = "\n".join(lines)
assert "calc_m_split" not in src, "strip failed"
assert "ellPQ" not in src, "strip failed"
open(PATH, "w").write(src)

# --- rebuild generator pieces in long double ---------------------------------
# reuse the section templates from gen_calc_m_split.py WITHOUT executing its
# splice logic: exec only up to the splice marker
class _G:
    pass


g = _G()
gen_src = open(os.path.join(CODEGEN, "gen_calc_m_split.py")).read()
gen_src = gen_src.split("# --- splice into the file")[0]
gen_src = gen_src.replace('sys.exit("calc_m_split already present; aborting")',
                          "pass")
exec(gen_src, g.__dict__)

series26 = open(os.path.join(CODEGEN, "pq_series_26.txt")).read()
# -> long double arrays with L-suffixed literals
import re
series26 = series26.replace("static const double", "static const long double")
series26 = re.sub(r"(-?\d+\.\d+(?:e-?\d+)?)", r"\1L", series26)
series26 = series26.replace("1.0L,", "1.0L,").replace("{\n  0.5L,", "{\n  0.5L,")

SERIES = ("\n/* Taylor coefficients in m1 = C1/(1+C1) for the log-split pieces\n"
          "   K(1-m1) = PK + QK*ln(1/m1), E(1-m1) = PE + QE*ln(1/m1), used for\n"
          "   the |m1| < 0.1 branch of the PQ primitive where the closed-form\n"
          "   Legendre-derivative expressions suffer 1/m1, 1/m1^2 cancellation.\n"
          "   Long double: the split assembly cancels ~C1^3 near the particle,\n"
          "   so every input factor needs better-than-double relative accuracy. */\n"
          + series26 + "\n")

PRIMITIVE = r"""
/* PK, QK, PE, QE and their first/second derivatives with respect to C1,
   in long double.  All twelve outputs are analytic in C1 >= 0 (no poles at
   C1 = 0): the 1/m1 poles of the classical derivative identities live
   entirely in the ln(1/m1) channel bookkeeping, which calc_m_split handles
   via ln(alpha).  Layout:
   out = {PK,QK,PE,QE, dPK,dQK,dPE,dQE, d2PK,d2QK,d2PE,d2QE}. */
static void effsource_ellPQ_derivs_from_C1l(long double C1, long double *out)
{
  const long double s = 1.0L / (1.0L + C1);
  const long double m1 = C1 * s;
  long double f[4], fm[4], fmm[4];

  if (m1 < 0.1L)
  {
    const long double *ser[4] = {effsource_pk_ser, effsource_qk_ser,
                                 effsource_pe_ser, effsource_qe_ser};
    for (int i = 0; i < 4; i++)
    {
      long double v = 0.0L, d1 = 0.0L, d2 = 0.0L;
      for (int n = EFFSOURCE_PQ_NSER - 1; n >= 0; n--)
      {
        v = v * m1 + ser[i][n];
        if (n >= 1) d1 = d1 * m1 + n * ser[i][n];
        if (n >= 2) d2 = d2 * m1 + n * (n - 1) * ser[i][n];
      }
      f[i] = v; fm[i] = d1; fmm[i] = d2;
    }
  }
  else
  {
    /* Away from m1 = 0 the cancellation is bounded (amplification <= 1/m1
       = 10) and double-precision AGM values suffice. */
    double Kmd, Emd, Kcd, Ecd;
    effsource_ellKE_from_C1(1.0 / (double)C1, &Kmd, &Emd); /* parameter m1 */
    effsource_ellKE_from_C1((double)C1, &Kcd, &Ecd);       /* param 1 - m1 */
    const long double Km = Kmd, Em = Emd, Kc = Kcd, Ec = Ecd;
    const long double L = -logl(m1);
    const long double m1c = 1.0L - m1;
    const long double QK = Km / M_PI, QE = (Km - Em) / M_PI;
    const long double PK = Kc - QK * L, PE = Ec - QE * L;
    /* Legendre derivative identities in the parameter m1 */
    const long double Kmp = (Em - m1c * Km) / (2.0L * m1 * m1c);
    const long double Emp = (Em - Km) / (2.0L * m1);
    const long double Kcp = -(Ec - m1 * Kc) / (2.0L * m1 * m1c);
    const long double Ecp = -(Ec - Kc) / (2.0L * m1c);
    const long double D = 2.0L * m1 * m1c, Dp = 2.0L * (1.0L - 2.0L * m1);
    const long double Kmpp = ((Emp + Km - m1c * Kmp) - Kmp * Dp) / D;
    const long double Empp = -(Emp + Kmp) / (2.0L * m1);
    const long double Kcpp = (-(Ecp - Kc - m1 * Kcp) - Kcp * Dp) / D;
    const long double Ecpp = (Kcp + Ecp) / (2.0L * m1c);
    const long double QKp = Kmp / M_PI, QEp = (Kmp - Emp) / M_PI;
    const long double PKp = Kcp - QKp * L + QK / m1;
    const long double PEp = Ecp - QEp * L + QE / m1;
    const long double QKpp = Kmpp / M_PI, QEpp = (Kmpp - Empp) / M_PI;
    const long double PKpp = Kcpp - QKpp * L + 2.0L * QKp / m1 - QK / (m1 * m1);
    const long double PEpp = Ecpp - QEpp * L + 2.0L * QEp / m1 - QE / (m1 * m1);
    f[0] = PK; f[1] = QK; f[2] = PE; f[3] = QE;
    fm[0] = PKp; fm[1] = QKp; fm[2] = PEp; fm[3] = QEp;
    fmm[0] = PKpp; fmm[1] = QKpp; fmm[2] = PEpp; fmm[3] = QEpp;
  }

  /* chain rule m1 -> C1: dm1/dC1 = s^2, d2m1/dC1^2 = -2 s^3 */
  const long double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
  for (int i = 0; i < 4; i++)
  {
    out[i] = f[i];
    out[4 + i] = fm[i] * s2;
    out[8 + i] = fmm[i] * s4 - 2.0L * fm[i] * s3;
  }
}

/* double wrapper (also handy for external validation) */
void effsource_ellPQ_derivs_from_C1(double C1, double *out)
{
  long double buf[12];
  effsource_ellPQ_derivs_from_C1l(C1, buf);
  for (int i = 0; i < 12; i++)
    out[i] = (double)buf[i];
}
"""


def to_ld(text):
    """Type-transform a section to long double."""
    text = text.replace("double _Complex", "long double _Complex")
    text = text.replace("const double ", "const long double ")
    text = re.sub(r"(?<![a-zA-Z_])double ", "long double ", text)
    text = text.replace("long long double", "long double")
    for fn in ("cexp", "creal", "cimag"):
        text = re.sub(rf"(?<![a-zA-Z_]){fn}\(", fn + "l(", text)
    for fn in ("sqrt", "log", "sin", "cos", "pow"):
        text = re.sub(rf"(?<![a-zA-Z_]){fn}\(", fn + "l(", text)
    return text


# fresh landmark discovery on the stripped file
lines = open(PATH).read().split("\n")


def find1(pat, start=0):
    for i in range(start, len(lines)):
        if pat in lines[i]:
            return i + 1
    raise SystemExit(f"landmark not found: {pat}")


core_sig = find1("static void effsource_equatorial_ctx_calc_m_core")
prelude_a = find1("const double M = ctx->M", core_sig)
ell_comment = find1("/* Elliptic integrals */", core_sig)
prelude_b = ell_comment - 1
while not lines[prelude_b - 1].strip():
    prelude_b -= 1
den_a = find1("/* Denominator - there is a different", core_sig)
rot_comment = find1("/* m-modes for the rotated", den_a)
den_b = rot_comment - 1
while not lines[den_b - 1].strip():
    den_b -= 1
i = rot_comment
while lines[i - 1] != "}":
    i += 1
core_end = i
pq_end = find1("*PE = ellE - *QE * lninv;")
pq_close = pq_end
while lines[pq_close - 1] != "}":
    pq_close += 1

PRELUDE = to_ld("\n".join(lines[prelude_a - 1:prelude_b]))
DEN = to_ld("\n".join(lines[den_a - 1:den_b]))
BOX = to_ld(g.BOX_HELPER)
ELLIP = g.ELLIP_PQ.replace("double PQ[12];", "long double PQ[12];") \
                  .replace("effsource_ellPQ_derivs_from_C1(C1, PQ);",
                           "effsource_ellPQ_derivs_from_C1l(C1, PQ);")
ELLIP = to_ld(ELLIP).replace("effsource_ellPQ_derivs_from_C1l(",
                             "effsource_ellPQ_derivs_from_C1l(")
NUM = to_ld(g.num_loops())
ASSEMBLY = to_ld(g.ASSEMBLY)
# output stores: creall into double arrays -- add explicit casts for clarity
ASSEMBLY = ASSEMBLY.replace("creall(", "(double)creall(") \
                   .replace("cimagl(", "(double)cimagl(")
# the output buffers (PhiS_s/dPhiS_s/d2PhiS_s/src_s) are double* -- keep the
# pointers into them double (to_ld wrongly promotes these declarations)
ASSEMBLY = ASSEMBLY.replace("long double *dP ", "double *dP ") \
                   .replace("long double *d2P ", "double *d2P ")
# the box helper returns ld complex; calls unaffected

HEADER = g.HEADER
ENTRY = g.ENTRY

func = (HEADER + PRELUDE + "\n" + ELLIP + "\n" + NUM + "\n"
        + DEN + "\n" + ASSEMBLY + ENTRY)

new_lines = list(lines)
new_lines[pq_close - 1] += "\n" + SERIES + PRIMITIVE
new_lines[core_end - 1] += "\n" + BOX + func
open(PATH, "w").write("\n".join(new_lines))
print(f"ld re-splice done: prelude {prelude_a}-{prelude_b}, den {den_a}-{den_b}")

import subprocess
subprocess.run(["gcc", "-fsyntax-only",
                "-I/home/jonathan/projects/sf/effectivesource", PATH],
               check=False)
