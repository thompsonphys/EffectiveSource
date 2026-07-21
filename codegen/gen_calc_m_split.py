"""Assemble effsource_equatorial_ctx_calc_m_split into kerr-equatorial-ctx.c.

Verbatim-copies the calc_m_core prelude (coefficient polynomials, alpha/C
power arrays) and denominator block so the giant generated expressions are
bit-identical, and splices in the P/Q elliptic split, doubled numerator
loops, and the complex-arithmetic channel assembly.

Idempotent: refuses to run if calc_m_split already present.
"""
import re
import subprocess
import sys

PATH = "/home/jonathan/projects/sf/effectivesource/kerr-equatorial-ctx.c"
src = open(PATH).read()
if "calc_m_split" in src:
    sys.exit("calc_m_split already present; aborting")

lines = src.split("\n")  # 0-based; file line N = lines[N-1]


def block(a, b):
    """File lines a..b inclusive (1-based)."""
    return "\n".join(lines[a - 1:b])


# --- landmark discovery (don't trust hardcoded numbers) -----------------------
def find(pat, start=0):
    for i in range(start, len(lines)):
        if pat in lines[i]:
            return i + 1
    raise SystemExit(f"landmark not found: {pat}")


core_sig = find("static void effsource_equatorial_ctx_calc_m_core")
prelude_a = find("const double M = ctx->M", core_sig)
ell_comment = find("/* Elliptic integrals */", core_sig)
prelude_b = ell_comment - 1
while not lines[prelude_b - 1].strip():
    prelude_b -= 1
den_a = find("/* Denominator - there is a different", core_sig)
rot_comment = find("/* m-modes for the rotated", den_a)
den_b = rot_comment - 1
while not lines[den_b - 1].strip():
    den_b -= 1
# closing brace of calc_m_core
i = rot_comment
while lines[i - 1] != "}":
    i += 1
core_end = i
pq_end = find("*PE = ellE - *QE * lninv;")
pq_close = pq_end
while lines[pq_close - 1] != "}":
    pq_close += 1

PRELUDE = block(prelude_a, prelude_b)
DEN = block(den_a, den_b)

# --- series coefficients (from gen_pq_series.py, mpmath dps=60) --------------
SERIES = r"""
/* Taylor coefficients in m1 = C1/(1+C1) for the log-split pieces
   K(1-m1) = PK + QK*ln(1/m1), E(1-m1) = PE + QE*ln(1/m1), used for the
   |m1| < 0.1 branch of effsource_ellPQ_derivs_from_C1 where the closed-form
   Legendre-derivative expressions suffer 1/m1, 1/m1^2 cancellation. */
#define EFFSOURCE_PQ_NSER 24
static const double effsource_pk_ser[EFFSOURCE_PQ_NSER] = {
  1.3862943611198906188, 0.096573590279972654709, 0.030885144532484618274,
  0.014937600369780984912, 0.0087663121971760665734, 0.0057548876844001139245,
  0.0040646585499939759365, 0.0030225465503440763496, 0.0023351571703480847014,
  0.0018580702734561832841, 0.0015135115737195792364, 0.0012565911382975057151,
  0.001059929706785451804, 0.00090605962449273640226, 0.00078341184444331734633,
  0.00068407788767568759869, 0.00060250319343517162842, 0.00053469388084383163545,
  0.00047771869488793861634, 0.00042938696968361928421, 0.00038803495047317503885,
  0.0003523807257523782792, 0.00032142370171640571553, 0.00029437364853540258112};
static const double effsource_qk_ser[EFFSOURCE_PQ_NSER] = {
  0.5, 0.125, 0.0703125, 0.048828125, 0.037384033203125, 0.03028106689453125,
  0.025444507598876953125, 0.021939396858215332031, 0.019282673019915819168,
  0.0171996682183817029, 0.015522700567089486867, 0.014143617665467900224,
  0.012989537751792568088, 0.012009557832648454223, 0.011167050586735616236,
  0.010434988381605170282, 0.0097929920260962584387, 0.0092254051180093645673,
  0.008720000979599900922, 0.0082670923414627869545, 0.0078589046571030618486,
  0.0074891262633731558773, 0.0071525797835624820336, 0.0068449782900349839877};
static const double effsource_pe_ser[EFFSOURCE_PQ_NSER] = {
  1.0, 0.44314718055994530942, 0.056805192709979491031, 0.021831370443737181895,
  0.011544521417308361798, 0.0071420003133959599161, 0.0048547433371649481808,
  0.0035146879637813760929, 0.0026622358529927642963, 0.0020863973706379082725,
  0.0016791684186914656054, 0.0013805722023228649452, 0.001155123390641123879,
  0.00098073259453692091186, 0.00084306372506532118834, 0.00073248244832468317649,
  0.00064231961700299976303, 0.00056783961849660038945, 0.00050560453674938976199,
  0.00045306958660196806976, 0.00040831843001470031841, 0.00036988570034047692358,
  0.00033663537220906683497, 0.00030767575755474237575};
static const double effsource_qe_ser[EFFSOURCE_PQ_NSER] = {
  0.0, 0.25, 0.09375, 0.05859375, 0.042724609375, 0.0336456298828125,
  0.0277576446533203125, 0.023627042770385742187, 0.020568184554576873779,
  0.018211413407698273659, 0.016339684807462617755, 0.014817123268585419282,
  0.013554300262740071048, 0.012489940145954392392, 0.011580645052911009429,
  0.010794815567177762361, 0.010108894994680008711, 0.0095049628488581331905,
  0.0089691438647313266626, 0.0084905272696104298451, 0.0080604150329262172806,
  0.0076717878795529889474, 0.0073189188482964932437, 0.0069970889187024280763};

/* PK, QK, PE, QE and their first/second derivatives with respect to C1.
   All twelve outputs are analytic in C1 >= 0 (no poles at C1 = 0): the
   1/m1 poles of the classical derivative identities live entirely in the
   ln(1/m1) channel bookkeeping, which the caller handles via ln(alpha).
   Layout: out = {PK,QK,PE,QE, dPK,dQK,dPE,dQE, d2PK,d2QK,d2PE,d2QE}. */
static void effsource_ellPQ_derivs_from_C1(double C1, double *out)
{
  const double s = 1.0 / (1.0 + C1);
  const double m1 = C1 * s;
  double f[4], fm[4], fmm[4];

  if (m1 < 0.1)
  {
    const double *ser[4] = {effsource_pk_ser, effsource_qk_ser,
                            effsource_pe_ser, effsource_qe_ser};
    for (int i = 0; i < 4; i++)
    {
      double v = 0.0, d1 = 0.0, d2 = 0.0;
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
    double Km, Em, Kc, Ec;
    effsource_ellKE_from_C1(1.0 / C1, &Km, &Em);   /* parameter m1     */
    effsource_ellKE_from_C1(C1, &Kc, &Ec);         /* parameter 1 - m1 */
    const double L = -log(m1);
    const double m1c = 1.0 - m1;
    const double QK = Km / M_PI, QE = (Km - Em) / M_PI;
    const double PK = Kc - QK * L, PE = Ec - QE * L;
    /* Legendre derivative identities in the parameter m1 */
    const double Kmp = (Em - m1c * Km) / (2.0 * m1 * m1c);
    const double Emp = (Em - Km) / (2.0 * m1);
    const double Kcp = -(Ec - m1 * Kc) / (2.0 * m1 * m1c);
    const double Ecp = -(Ec - Kc) / (2.0 * m1c);
    const double D = 2.0 * m1 * m1c, Dp = 2.0 * (1.0 - 2.0 * m1);
    const double Kmpp = ((Emp + Km - m1c * Kmp) - Kmp * Dp) / D;
    const double Empp = -(Emp + Kmp) / (2.0 * m1);
    const double Kcpp = (-(Ecp - Kc - m1 * Kcp) - Kcp * Dp) / D;
    const double Ecpp = (Kcp + Ecp) / (2.0 * m1c);
    const double QKp = Kmp / M_PI, QEp = (Kmp - Emp) / M_PI;
    const double PKp = Kcp - QKp * L + QK / m1;
    const double PEp = Ecp - QEp * L + QE / m1;
    const double QKpp = Kmpp / M_PI, QEpp = (Kmpp - Empp) / M_PI;
    const double PKpp = Kcpp - QKpp * L + 2.0 * QKp / m1 - QK / (m1 * m1);
    const double PEpp = Ecpp - QEpp * L + 2.0 * QEp / m1 - QE / (m1 * m1);
    f[0] = PK; f[1] = QK; f[2] = PE; f[3] = QE;
    fm[0] = PKp; fm[1] = QKp; fm[2] = PEp; fm[3] = QEp;
    fmm[0] = PKpp; fmm[1] = QKpp; fmm[2] = PEpp; fmm[3] = QEpp;
  }

  /* chain rule m1 -> C1: dm1/dC1 = s^2, d2m1/dC1^2 = -2 s^3 */
  const double s2 = s * s, s3 = s2 * s, s4 = s2 * s2;
  for (int i = 0; i < 4; i++)
  {
    out[i] = f[i];
    out[4 + i] = fm[i] * s2;
    out[8 + i] = fmm[i] * s4 - 2.0 * fm[i] * s3;
  }
}
"""

BOX_HELPER = r"""
/* Scalar wave operator on one channel of the m-mode field (complex-valued
   copy of the Box[PhiS] expression in calc_m_core). */
static double _Complex effsource_box_m_channel(double a, double r, double theta,
    double _Complex dPr, double _Complex dPth,
    double _Complex dtt, double _Complex dtph, double _Complex drr,
    double _Complex dthth, double _Complex dphph)
{
  const double sinth = sin(theta);
  const double sinth2 = sinth * sinth;
  const double sin2th = sin(2.0 * theta);
  const double cos2th = cos(2.0 * theta);
  const double r2 = r * r;
  const double r3 = r2 * r;
  const double r4 = r2 * r2;
  const double a2 = a * a;
  const double a4 = a2 * a2;

  return -((2 * a2 * dPr - a2 * dphph - a4 * drr - a2 * dthth - 4 * dPr * r - 2 * a2 * dPr * r +
            4 * dphph * r + 2 * a * dtph * r + 4 * a2 * drr * r + 2 * dthth * r + 6 * dPr * r2 - 2 * dphph * r2 - 4 * drr * r2 -
            2 * a2 * drr * r2 - dthth * r2 - 2 * dPr * r3 + 4 * drr * r3 - drr * r4 +
            (a4 * dtt + 4 * a * dtph * r + 2 * dtt * r4 + a2 * dtt * r * (2 + 3 * r)) * sinth2 +
            cos2th * (a4 * drr - 2 * a * dtph * r + (-2 + r) * r * (dthth + 2 * dPr * (-1 + r) + drr * (-2 + r) * r) +
                      a2 * (-dphph + dthth + 2 * dPr * (-1 + r) - 4 * drr * r + 2 * drr * r2) +
                      a2 * dtt * (a2 + (-2 + r) * r) * sinth2) -
            a2 * dPth * sin2th + 2 * dPth * r * sin2th -
            dPth * r2 * sin2th)) /
         ((sinth2 * (a2 + (-2 + r) * r) * (a2 + 2 * r2 + a2 * cos2th)));
}
"""

ELLIP_PQ = r"""
  /* Elliptic integrals: log-split channels.  K = PK + QK*ln(1/m1) and
     E = PE + QE*ln(1/m1) with ln(1/m1) = ln(alpha+beta) - ln(alpha); the
     P/Q pieces and their C1-derivatives are analytic through the particle,
     so every 1/alpha, 1/alpha^2 pole of the classical dK/dC identities is
     deferred to the explicit ln(alpha) bookkeeping below. */
  double PQ[12];
  effsource_ellPQ_derivs_from_C1(C1, PQ);
  const double ellipP[2] = {PQ[0], PQ[2]};
  const double ellipQ[2] = {PQ[1], PQ[3]};
  const double dellipP_dC[2] = {PQ[4], PQ[6]};
  const double dellipQ_dC[2] = {PQ[5], PQ[7]};
  const double d2ellipP_dC2[2] = {PQ[8], PQ[10]};
  const double d2ellipQ_dC2[2] = {PQ[9], PQ[11]};

  double dellipP_dr[2], dellipP_dtheta[2], dellipP_dt[2];
  double d2ellipP_dr2[2], d2ellipP_dtheta2[2], d2ellipP_dt2[2];
  double dellipQ_dr[2], dellipQ_dtheta[2], dellipQ_dt[2];
  double d2ellipQ_dr2[2], d2ellipQ_dtheta2[2], d2ellipQ_dt2[2];
  for (int i = 0; i < 2; i++)
  {
    dellipP_dr[i] = dellipP_dC[i] * dC1_dr;
    dellipP_dtheta[i] = dellipP_dC[i] * dC1_dtheta;
    dellipP_dt[i] = dellipP_dC[i] * dC1_dt;
    d2ellipP_dr2[i] = d2ellipP_dC2[i] * dC1_dr * dC1_dr + dellipP_dC[i] * d2C1_dr2;
    d2ellipP_dtheta2[i] = d2ellipP_dC2[i] * dC1_dtheta * dC1_dtheta + dellipP_dC[i] * d2C1_dtheta2;
    d2ellipP_dt2[i] = d2ellipP_dC2[i] * dC1_dt * dC1_dt + dellipP_dC[i] * d2C1_dt2;
    dellipQ_dr[i] = dellipQ_dC[i] * dC1_dr;
    dellipQ_dtheta[i] = dellipQ_dC[i] * dC1_dtheta;
    dellipQ_dt[i] = dellipQ_dC[i] * dC1_dt;
    d2ellipQ_dr2[i] = d2ellipQ_dC2[i] * dC1_dr * dC1_dr + dellipQ_dC[i] * d2C1_dr2;
    d2ellipQ_dtheta2[i] = d2ellipQ_dC2[i] * dC1_dtheta * dC1_dtheta + dellipQ_dC[i] * d2C1_dtheta2;
    d2ellipQ_dt2[i] = d2ellipQ_dC2[i] * dC1_dt * dC1_dt + dellipQ_dC[i] * d2C1_dt2;
  }

  if (m > 20)
  {
    printf("Support for computing mode %d has not yet been added.\n", m);
    return;
  }
"""


def num_loops():
    """Emit the seven numerator loop nests with P and Q accumulators."""
    out = []

    def loop(suffix, facsP_re, facsQ_re, facsP_im, facsQ_im, weights):
        """weights: list of C-array factors pairing with the fac lists."""
        vP_re, vP_im = f"NumP_re{suffix}", f"NumP_im{suffix}"
        vQ_re, vQ_im = f"NumQ_re{suffix}", f"NumQ_im{suffix}"
        s = [f"  double {vP_re} = 0, {vP_im} = 0, {vQ_re} = 0, {vQ_im} = 0;"]
        s.append("  for (int i = 0; i < 2; i++)")
        s.append("    for (int j = 0; j < 5; j++)")
        s.append("    {")
        for n, (fp, fq) in enumerate(zip(facsP_re, facsQ_re)):
            s.append(f"      const double fPre{n} = {fp};")
            s.append(f"      const double fQre{n} = {fq};")
        s.append("      for (int k = max(j - i, 0); k <= m + 2 + j; k++)")
        s.append("      {")
        s.append("        const double e = ReEI[m][i][j][k];")
        termP = " + ".join(f"{w} * fPre{n}" for n, w in enumerate(weights))
        termQ = " + ".join(f"{w} * fQre{n}" for n, w in enumerate(weights))
        s.append(f"        {vP_re} += e * ({termP});")
        s.append(f"        {vQ_re} += e * ({termQ});")
        s.append("      }")
        for n, (fp, fq) in enumerate(zip(facsP_im, facsQ_im)):
            s.append(f"      const double fPim{n} = {fp};")
            s.append(f"      const double fQim{n} = {fq};")
        s.append("      for (int k = max(j - i - 1, 0); k <= m + 1 + j; k++)")
        s.append("      {")
        s.append("        const double e = ImEI[m][i][j][k];")
        termP = " + ".join(f"{w} * fPim{n}" for n, w in enumerate(weights))
        termQ = " + ".join(f"{w} * fQim{n}" for n, w in enumerate(weights))
        s.append(f"        {vP_im} += e * ({termP});")
        s.append(f"        {vQ_im} += e * ({termQ});")
        s.append("      }")
        s.append("    }")
        s.append("")
        return "\n".join(s)

    # value
    out.append(loop("",
        ["ellipP[i] * ReA[j]"], ["ellipQ[i] * ReA[j]"],
        ["ellipP[i] * ImA[j]"], ["ellipQ[i] * ImA[j]"],
        ["C[k]"]))
    # first derivatives
    for d in ("r", "theta", "t"):
        out.append(loop(f"_d{d}",
            [f"dellipP_d{d}[i] * ReA[j] + ellipP[i] * dReA_d{d}[j]",
             "ellipP[i] * ReA[j]"],
            [f"dellipQ_d{d}[i] * ReA[j] + ellipQ[i] * dReA_d{d}[j]",
             "ellipQ[i] * ReA[j]"],
            [f"dellipP_d{d}[i] * ImA[j] + ellipP[i] * dImA_d{d}[j]",
             "ellipP[i] * ImA[j]"],
            [f"dellipQ_d{d}[i] * ImA[j] + ellipQ[i] * dImA_d{d}[j]",
             "ellipQ[i] * ImA[j]"],
            ["C[k]", f"dC_d{d}[k]"]))
    # second derivatives
    for d in ("r", "theta", "t"):
        out.append(loop(f"_d{d}2",
            [f"d2ellipP_d{d}2[i] * ReA[j] + 2.0 * dellipP_d{d}[i] * dReA_d{d}[j] + ellipP[i] * d2ReA_d{d}2[j]",
             f"2.0 * (dellipP_d{d}[i] * ReA[j] + ellipP[i] * dReA_d{d}[j])",
             "ellipP[i] * ReA[j]"],
            [f"d2ellipQ_d{d}2[i] * ReA[j] + 2.0 * dellipQ_d{d}[i] * dReA_d{d}[j] + ellipQ[i] * d2ReA_d{d}2[j]",
             f"2.0 * (dellipQ_d{d}[i] * ReA[j] + ellipQ[i] * dReA_d{d}[j])",
             "ellipQ[i] * ReA[j]"],
            [f"d2ellipP_d{d}2[i] * ImA[j] + 2.0 * dellipP_d{d}[i] * dImA_d{d}[j] + ellipP[i] * d2ImA_d{d}2[j]",
             f"2.0 * (dellipP_d{d}[i] * ImA[j] + ellipP[i] * dImA_d{d}[j])",
             "ellipP[i] * ImA[j]"],
            [f"d2ellipQ_d{d}2[i] * ImA[j] + 2.0 * dellipQ_d{d}[i] * dImA_d{d}[j] + ellipQ[i] * d2ImA_d{d}2[j]",
             f"2.0 * (dellipQ_d{d}[i] * ImA[j] + ellipQ[i] * dImA_d{d}[j])",
             "ellipQ[i] * ImA[j]"],
            ["C[k]", f"dC_d{d}[k]", f"d2C_d{d}2[k]"]))
    return "\n".join(out)


# NB: original loops name the theta/t derivative arrays dReA_dtheta / dReA_dt
# and dC_dtheta / dC_dt, and second derivatives d2ReA_dr2 / d2ReA_dtheta2 /
# d2ReA_dt2, d2C_dr2 / d2C_dtheta2 / d2C_dt2 -- the f-string patterns above
# reproduce exactly these names for d in {r, theta, t}.

ASSEMBLY = r"""
  /* Quotient rule per channel and part: f = N/D, f_x = (N_x - f D_x)/D,
     f_xx = (N_xx - 2 f_x D_x - f D_xx)/D.  Real parts use DenRePhiSb,
     imaginary parts DenImPhiSb, as in calc_m_core. */
  const double invDre = 1.0 / DenRePhiSb;
  const double invDim = 1.0 / DenImPhiSb;

  const double fPre = NumP_re * invDre;
  const double fPre_dr = (NumP_re_dr - fPre * dDenRePhiSb_dr) * invDre;
  const double fPre_dth = (NumP_re_dtheta - fPre * dDenRePhiSb_dtheta) * invDre;
  const double fPre_dt = (NumP_re_dt - fPre * dDenRePhiSb_dt) * invDre;
  const double fPre_dr2 = (NumP_re_dr2 - 2.0 * fPre_dr * dDenRePhiSb_dr - fPre * d2DenRePhiSb_dr2) * invDre;
  const double fPre_dth2 = (NumP_re_dtheta2 - 2.0 * fPre_dth * dDenRePhiSb_dtheta - fPre * d2DenRePhiSb_dtheta2) * invDre;
  const double fPre_dt2 = (NumP_re_dt2 - 2.0 * fPre_dt * dDenRePhiSb_dt - fPre * d2DenRePhiSb_dt2) * invDre;

  const double fPim = NumP_im * invDim;
  const double fPim_dr = (NumP_im_dr - fPim * dDenImPhiSb_dr) * invDim;
  const double fPim_dth = (NumP_im_dtheta - fPim * dDenImPhiSb_dtheta) * invDim;
  const double fPim_dt = (NumP_im_dt - fPim * dDenImPhiSb_dt) * invDim;
  const double fPim_dr2 = (NumP_im_dr2 - 2.0 * fPim_dr * dDenImPhiSb_dr - fPim * d2DenImPhiSb_dr2) * invDim;
  const double fPim_dth2 = (NumP_im_dtheta2 - 2.0 * fPim_dth * dDenImPhiSb_dtheta - fPim * d2DenImPhiSb_dtheta2) * invDim;
  const double fPim_dt2 = (NumP_im_dt2 - 2.0 * fPim_dt * dDenImPhiSb_dt - fPim * d2DenImPhiSb_dt2) * invDim;

  const double fQre = NumQ_re * invDre;
  const double fQre_dr = (NumQ_re_dr - fQre * dDenRePhiSb_dr) * invDre;
  const double fQre_dth = (NumQ_re_dtheta - fQre * dDenRePhiSb_dtheta) * invDre;
  const double fQre_dt = (NumQ_re_dt - fQre * dDenRePhiSb_dt) * invDre;
  const double fQre_dr2 = (NumQ_re_dr2 - 2.0 * fQre_dr * dDenRePhiSb_dr - fQre * d2DenRePhiSb_dr2) * invDre;
  const double fQre_dth2 = (NumQ_re_dtheta2 - 2.0 * fQre_dth * dDenRePhiSb_dtheta - fQre * d2DenRePhiSb_dtheta2) * invDre;
  const double fQre_dt2 = (NumQ_re_dt2 - 2.0 * fQre_dt * dDenRePhiSb_dt - fQre * d2DenRePhiSb_dt2) * invDre;

  const double fQim = NumQ_im * invDim;
  const double fQim_dr = (NumQ_im_dr - fQim * dDenImPhiSb_dr) * invDim;
  const double fQim_dth = (NumQ_im_dtheta - fQim * dDenImPhiSb_dtheta) * invDim;
  const double fQim_dt = (NumQ_im_dt - fQim * dDenImPhiSb_dt) * invDim;
  const double fQim_dr2 = (NumQ_im_dr2 - 2.0 * fQim_dr * dDenImPhiSb_dr - fQim * d2DenImPhiSb_dr2) * invDim;
  const double fQim_dth2 = (NumQ_im_dtheta2 - 2.0 * fQim_dth * dDenImPhiSb_dtheta - fQim * d2DenImPhiSb_dtheta2) * invDim;
  const double fQim_dt2 = (NumQ_im_dt2 - 2.0 * fQim_dt * dDenImPhiSb_dt - fQim * d2DenImPhiSb_dt2) * invDim;

  /* Complex rotated-frame quantities u = fRe + i fIm and the combined
     rotation+phase factor V = exp(-i m (c dr + phi_p)); psi collects the
     phase argument so V-derivatives are just products. */
  const double _Complex uP = fPre + I * fPim;
  const double _Complex uP_r = fPre_dr + I * fPim_dr;
  const double _Complex uP_th = fPre_dth + I * fPim_dth;
  const double _Complex uP_t = fPre_dt + I * fPim_dt;
  const double _Complex uP_rr = fPre_dr2 + I * fPim_dr2;
  const double _Complex uP_thth = fPre_dth2 + I * fPim_dth2;
  const double _Complex uP_tt = fPre_dt2 + I * fPim_dt2;

  const double _Complex uQ = fQre + I * fQim;
  const double _Complex uQ_r = fQre_dr + I * fQim_dr;
  const double _Complex uQ_th = fQre_dth + I * fQim_dth;
  const double _Complex uQ_t = fQre_dt + I * fQim_dt;
  const double _Complex uQ_rr = fQre_dr2 + I * fQim_dr2;
  const double _Complex uQ_thth = fQre_dth2 + I * fQim_dth2;
  const double _Complex uQ_tt = fQre_dt2 + I * fQim_dt2;

  const double psi_r = c;
  const double psi_t = dcdt * dr - c * rt + phit;
  const double psi_tt = d2cdt2 * dr - 2.0 * dcdt * rt - c * rtt + phitt;
  const double _Complex V = cexp(-I * (double)m * (c * dr + xp.phi));
  const double _Complex imr = I * (double)m * psi_r;
  const double _Complex imt = I * (double)m * psi_t;
  const double _Complex imtt = I * (double)m * psi_tt;
  const double m2r = (double)m * (double)m * psi_r * psi_r;
  const double m2t = (double)m * (double)m * psi_t * psi_t;

  const double _Complex GP = V * uP;
  const double _Complex GP_r = V * (uP_r - imr * uP);
  const double _Complex GP_th = V * uP_th;
  const double _Complex GP_t = V * (uP_t - imt * uP);
  const double _Complex GP_rr = V * (uP_rr - 2.0 * imr * uP_r - m2r * uP);
  const double _Complex GP_thth = V * uP_thth;
  const double _Complex GP_tt = V * (uP_tt - 2.0 * imt * uP_t - (imtt + m2t) * uP);

  const double _Complex GQ = V * uQ;
  const double _Complex GQ_r = V * (uQ_r - imr * uQ);
  const double _Complex GQ_th = V * uQ_th;
  const double _Complex GQ_t = V * (uQ_t - imt * uQ);
  const double _Complex GQ_rr = V * (uQ_rr - 2.0 * imr * uQ_r - m2r * uQ);
  const double _Complex GQ_thth = V * uQ_thth;
  const double _Complex GQ_tt = V * (uQ_tt - 2.0 * imt * uQ_t - (imtt + m2t) * uQ);

  /* ln(1/m1) = lnR - ln(alpha) with R = alpha + beta.  lnR is analytic. */
  const double Rab = alpha_plus_beta_10;
  const double lnR = log(Rab);
  const double lnR_r = dalpha_dr / Rab;
  const double lnR_th = dalpha_dtheta / Rab;
  const double lnR_t = (dalpha_dt + dbetadt) / Rab;
  const double lnR_rr = d2alpha_dr2 / Rab - lnR_r * lnR_r;
  const double lnR_thth = d2alpha_dtheta2 / Rab - lnR_th * lnR_th;
  const double lnR_tt = (d2alpha_dt2 + d2betadt2) / Rab - lnR_t * lnR_t;

  /* Channels: value = A + L*ln(alpha) + P1/alpha + P2/alpha^2 */
  const double _Complex vA = GP + GQ * lnR;
  const double _Complex vL = -GQ;

  const double _Complex rA = GP_r + GQ_r * lnR + GQ * lnR_r;
  const double _Complex rL = -GQ_r;
  const double _Complex rP1 = -GQ * dalpha_dr;

  const double _Complex thA = GP_th + GQ_th * lnR + GQ * lnR_th;
  const double _Complex thL = -GQ_th;
  const double _Complex thP1 = -GQ * dalpha_dtheta;

  const double _Complex tA = GP_t + GQ_t * lnR + GQ * lnR_t;
  const double _Complex tL = -GQ_t;
  const double _Complex tP1 = -GQ * dalpha_dt;

  const double _Complex rrA = GP_rr + GQ_rr * lnR + 2.0 * GQ_r * lnR_r + GQ * lnR_rr;
  const double _Complex rrL = -GQ_rr;
  const double _Complex rrP1 = -(2.0 * GQ_r * dalpha_dr + GQ * d2alpha_dr2);
  const double _Complex rrP2 = GQ * dalpha_dr * dalpha_dr;

  const double _Complex ththA = GP_thth + GQ_thth * lnR + 2.0 * GQ_th * lnR_th + GQ * lnR_thth;
  const double _Complex ththL = -GQ_thth;
  const double _Complex ththP1 = -(2.0 * GQ_th * dalpha_dtheta + GQ * d2alpha_dtheta2);
  const double _Complex ththP2 = GQ * dalpha_dtheta * dalpha_dtheta;

  const double _Complex ttA = GP_tt + GQ_tt * lnR + 2.0 * GQ_t * lnR_t + GQ * lnR_tt;
  const double _Complex ttL = -GQ_tt;
  const double _Complex ttP1 = -(2.0 * GQ_t * dalpha_dt + GQ * d2alpha_dt2);
  const double _Complex ttP2 = GQ * dalpha_dt * dalpha_dt;

  /* phi dependence: d/dphi = i m, so channels scale */
  const double _Complex im = I * (double)m;
  const double mm = (double)m * (double)m;

  /* --- P-side hidden-pole closed forms (analytic in dr, dtheta, worldline
     data; generated, see codegen/poles_calc_m.txt).  These are the alpha^-1..
     alpha^-5 Laurent coefficients of the P-side ratio uP for each component;
     multiplied by V and moved out of the A channel into P1..P5 so A is analytic.
     Recombination A + L*lnalpha + sum_q Pq/alpha^q stays exact by construction. */
__POLE_DEFS__

  const double ia1 = 1.0 / alpha, ia2 = ia1 * ia1, ia3 = ia2 * ia1,
               ia4 = ia2 * ia2, ia5 = ia3 * ia2;

#define PS5(c) {V*(P1re_##c + I*P1im_##c), V*(P2re_##c + I*P2im_##c), \
                V*(P3re_##c + I*P3im_##c), V*(P4re_##c + I*P4im_##c), \
                V*(P5re_##c + I*P5im_##c)}
  const double _Complex Pv[5] = PS5(val), Pr[5] = PS5(dr), Pth[5] = PS5(dth),
      Pt[5] = PS5(dt), Prr[5] = PS5(drr), Pthth[5] = PS5(dthth), Ptt[5] = PS5(dtt);
#undef PS5
#define PSUM(P) (P[0]*ia1 + P[1]*ia2 + P[2]*ia3 + P[3]*ia4 + P[4]*ia5)

  /* Seven channels {A, L, P1, P2, P3, P4, P5}.  A is cleaned of the P-side
     poles; P1/P2 add the P-side to the existing GQ-side (d ln alpha) content;
     P3/P4/P5 are P-side only. */
  const double _Complex vCh[7] = {vA - PSUM(Pv), vL, Pv[0], Pv[1], Pv[2], Pv[3], Pv[4]};
  const double _Complex tCh[7] = {tA - PSUM(Pt), tL, tP1 + Pt[0], Pt[1], Pt[2], Pt[3], Pt[4]};
  const double _Complex rCh[7] = {rA - PSUM(Pr), rL, rP1 + Pr[0], Pr[1], Pr[2], Pr[3], Pr[4]};
  const double _Complex thCh[7] = {thA - PSUM(Pth), thL, thP1 + Pth[0], Pth[1], Pth[2], Pth[3], Pth[4]};
  const double _Complex ttCh[7] = {ttA - PSUM(Ptt), ttL, ttP1 + Ptt[0], ttP2 + Ptt[1], Ptt[2], Ptt[3], Ptt[4]};
  const double _Complex rrCh[7] = {rrA - PSUM(Prr), rrL, rrP1 + Prr[0], rrP2 + Prr[1], Prr[2], Prr[3], Prr[4]};
  const double _Complex ththCh[7] = {ththA - PSUM(Pthth), ththL, ththP1 + Pthth[0], ththP2 + Pthth[1], Pthth[2], Pthth[3], Pthth[4]};
#undef PSUM

  /* Box[PhiS] per channel (channel-wise; Box is linear and ln alpha, alpha^-q
     are common scalars, so the split is recombination-exact for src too). */
  double _Complex boxCh[7];
  for (int ch = 0; ch < 7; ch++)
    boxCh[ch] = effsource_box_m_channel(a, r, theta, rCh[ch], thCh[ch],
        ttCh[ch], im * tCh[ch], rrCh[ch], ththCh[ch], -mm * vCh[ch]);

  for (int ch = 0; ch < 7; ch++)
  {
    PhiS_s[2 * ch] = creal(vCh[ch]);
    PhiS_s[2 * ch + 1] = cimag(vCh[ch]);

    double *dP = dPhiS_s + 8 * ch;
    dP[0] = creal(tCh[ch]);
    dP[1] = cimag(tCh[ch]);
    dP[2] = creal(rCh[ch]);
    dP[3] = cimag(rCh[ch]);
    dP[4] = creal(thCh[ch]);
    dP[5] = cimag(thCh[ch]);
    dP[6] = creal(im * vCh[ch]);
    dP[7] = cimag(im * vCh[ch]);

    double *d2P = d2PhiS_s + 20 * ch;
    for (int q = 0; q < 20; q++)
      d2P[q] = 0.0;
    d2P[0] = creal(ttCh[ch]);
    d2P[1] = cimag(ttCh[ch]);
    d2P[6] = creal(im * tCh[ch]);
    d2P[7] = cimag(im * tCh[ch]);
    d2P[8] = creal(rrCh[ch]);
    d2P[9] = cimag(rrCh[ch]);
    d2P[14] = creal(ththCh[ch]);
    d2P[15] = cimag(ththCh[ch]);
    d2P[18] = creal(-mm * vCh[ch]);
    d2P[19] = cimag(-mm * vCh[ch]);

    src_s[2 * ch] = creal(boxCh[ch]);
    src_s[2 * ch + 1] = cimag(boxCh[ch]);
  }
}
"""

# --- inject the generated P-side pole closed forms into the ASSEMBLY ----------
def _load_pole_defs():
    path = "/home/jonathan/projects/sf/effectivesource/codegen/poles_calc_m.txt"
    have = {}
    for ln in open(path):
        ln = ln.strip()
        if not ln or "=" not in ln:
            continue
        lhs, rhs = ln.split("=", 1)
        have[lhs.strip()] = rhs.rstrip(";").strip()
    comps = ("val", "dr", "dth", "dt", "drr", "dthth", "dtt")
    out = []
    for c in comps:
        for q in (1, 2, 3, 4, 5):
            for part in ("re", "im"):
                name = f"P{q}{part}_{c}"
                out.append(f"  const double {name} = {have.get(name, '0')};")
    return "\n".join(out)


ASSEMBLY = ASSEMBLY.replace("__POLE_DEFS__", _load_pole_defs())

HEADER = r"""
/* Kernel-split version of calc_m: every output component is returned as
   seven channels {A, L, P1, P2, P3, P4, P5} with

     value = A + L*ln(alpha) + P1/alpha + P2/alpha^2 + ... + P5/alpha^5,
     alpha = alpha20*dr^2 + alpha02*dtheta^2,

   and all channel coefficients analytic in (dr, dtheta) across the particle
   (the P-side hidden poles are extracted in closed form into P1..P5 so A is
   analytic; the field reaches 1/alpha^3, 1st derivatives 1/alpha^4, 2nd
   derivatives / src 1/alpha^5).  Output layout (Re/Im interleaved as in calc_m):
     PhiS_s[14]     : 7 channels x (Re, Im)
     dPhiS_s[56]    : 7 channels x dPhiS_dx[8]
     d2PhiS_s[140]  : 7 channels x d2PhiS_dx2[20]   (unavailable mixed
                      derivatives stored as 0 instead of calc_m's NAN)
     src_s[14]      : 7 channels x (Re, Im)
   See docs/nmode-kernel-factorization.md (Stage 2, sec 6). */
static void effsource_equatorial_ctx_calc_m_split_core(struct effsource_equatorial_ctx *ctx,
    int m, double r, double theta, double dr, double dtheta,
    double *PhiS_s, double *dPhiS_s, double *d2PhiS_s, double *src_s)
{
"""

ENTRY = r"""
void effsource_equatorial_ctx_calc_m_split(struct effsource_equatorial_ctx *ctx,
    int m, double dr, double dtheta,
    double *PhiS_s, double *dPhiS_s, double *d2PhiS_s, double *src_s)
{
  const struct coordinate xp = ctx->xp;
  effsource_equatorial_ctx_calc_m_split_core(ctx, m, xp.r + dr, xp.theta + dtheta,
      dr, dtheta, PhiS_s, dPhiS_s, d2PhiS_s, src_s);
}
"""

func = (HEADER + PRELUDE + "\n" + ELLIP_PQ + "\n" + num_loops() + "\n"
        + DEN + "\n" + ASSEMBLY + ENTRY)

# --- splice into the file ------------------------------------------------
new_lines = list(lines)
# 1) series + primitive after effsource_ellKE_PQ_from_C1's closing brace
new_lines[pq_close - 1] += "\n" + SERIES
# 2) box helper + split function after calc_m_core's closing brace
new_lines[core_end - 1] += "\n" + BOX_HELPER + func
# 3) complex.h include
for i, ln in enumerate(new_lines):
    if ln.startswith("#include <math.h>"):
        new_lines[i] = "#include <math.h>\n#include <complex.h>"
        break

open(PATH, "w").write("\n".join(new_lines))
print(f"spliced: prelude {prelude_a}-{prelude_b}, den {den_a}-{den_b}, "
      f"core_end {core_end}, pq_close {pq_close}")
subprocess.run(["gcc", "-fsyntax-only", "-I/home/jonathan/projects/sf/effectivesource",
                PATH], check=False)
