#!/usr/bin/env python3
"""Generate the four transpiled kernels into src/_generated_kernels.jl.

Each function's C body is extracted between its opening `{` and matching `}`,
transpiled, and wrapped with a Julia signature, output-container preamble, and a
return that assembles the results (Complex for the m-mode kernels)."""
from pathlib import Path
from c2jl import _match_paren
from transpile import transpile_body

HERE = Path(__file__).resolve().parent
SRC = (HERE.parent.parent / "kerr-circular-ctx.c").read_text()
OUT = HERE.parent / "src" / "_generated_kernels.jl"


def body_after(sig_marker):
    """Extract the { ... } body following the signature containing sig_marker."""
    start = SRC.index(sig_marker)
    brace = SRC.index("{", start)
    close = _match_paren_brace(SRC, brace)
    return SRC[brace + 1 : close]


def _match_paren_brace(s, open_idx):
    depth = 0
    for i in range(open_idx, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced braces")


def indent(lines, pad="    "):
    return "\n".join(pad + l if l else l for l in lines)


# ---------------------------------------------------------------------------
# Coordinate-preamble handling: the C body derives dr = r - rp etc. The offset
# variant of each kernel replaces those catastrophic subtractions with the
# offsets supplied directly, reconstructing only the benign absolute r, theta
# needed by the box operator.
# ---------------------------------------------------------------------------

_COORD_LHS = {"r", "theta", "phi", "rp", "thetap", "phip", "dr", "dtheta", "dphi"}
_COORD_RHS = {"x.r", "x.theta", "x.phi", "xp.r", "xp.theta", "xp.phi",
              "r - rp", "theta - thetap", "phi - phip"}


def _is_coord_line(line):
    s = line.strip()
    if "=" not in s:
        return False
    lhs, rhs = s.split("=", 1)
    return lhs.strip() in _COORD_LHS and rhs.strip() in _COORD_RHS


def offset_body(body):
    """Return (core, has_dphi): body with coordinate-preamble lines removed and
    the offset preamble spliced in at their original position."""
    idxs = [i for i, l in enumerate(body) if _is_coord_line(l)]
    first = idxs[0]
    idxset = set(idxs)
    has_dphi = any(body[i].strip().startswith("dphi ") for i in idxs)
    core = [l for i, l in enumerate(body) if i not in idxset]
    pre = ["rp = xp.r", "thetap = xp.theta", "dr = off.dr", "dtheta = off.dtheta"]
    if has_dphi:
        pre.append("dphi = off.dphi")
    pre += ["r = rp + dr", "theta = thetap + dtheta"]
    # lines 0..first-1 are all retained (first is the smallest coord index)
    return core[:first] + pre + core[first:], has_dphi


def build(name, doc, extra_args, extra_call, setup, retblock, sig_marker, comp):
    """Emit the absolute (Coordinate) and offset (Offset) methods for a kernel,
    sharing an identical computational core.

    If `comp`, the kernel contains an opt-in compensated sum: emit a `_name`
    core taking a `cm::Val` mode plus a public `name(...; comp::Bool=false)`
    keyword wrapper. Otherwise emit the public method directly (no compensation).
    """
    body = transpile_body(body_after(sig_marker))
    off, _ = offset_body(body)
    setup_block = ("\n" + indent(setup)) if setup else ""
    ea = extra_args   # signature extras, e.g. "m::Integer, "
    eac = extra_call  # call-site forwarding, e.g. "m, "

    def core(argname, argtype, b):
        fname = f"_{name}" if comp else name
        cmarg = ", cm::Val" if comp else ""
        return (
            f"function {fname}(ctx::EffsourceCtx{{T}}, {ea}{argname}::{argtype}{{T}}{cmarg}) where {{T}}"
            + setup_block + "\n" + indent(b) + "\n" + retblock + "\nend\n"
        )

    def wrapper(argname, argtype):
        return (
            f"{name}(ctx::EffsourceCtx{{T}}, {ea}{argname}::{argtype}{{T}}; comp::Bool=false) where {{T}} = "
            f"_{name}(ctx, {eac}{argname}, Val(comp))\n"
        )

    offdoc = f'"""{name} at an `Offset` from the particle (avoids the r-rp cancellation)."""\n'
    out = [f'"""{doc}"""\n']
    if comp:
        out += [wrapper("x", "Coordinate"), core("x", "Coordinate", body),
                "\n", offdoc, wrapper("off", "Offset"), core("off", "Offset", off)]
    else:
        out += [core("x", "Coordinate", body), "\n", offdoc, core("off", "Offset", off)]
    return "".join(out)


def gen_phis():
    doc = "Real-space singular field PhiS at field point `x` (or `Offset`)."
    return build("phis", doc, "", "", [], "    return PhiS_out",
                 "void effsource_ctx_PhiS(struct", comp=False)


def gen_calc():
    doc = ("Singular field, its first/second derivatives and the effective source.\n"
           "Returns `(PhiS, dPhiS_dx[4], d2PhiS_dx2[10], src)` with component order\n"
           "(t,r,theta,phi) and (tt,tr,ttheta,tphi,rr,rtheta,rphi,thetatheta,thetaphi,phiphi).\n"
           "`comp=true` sums the source numerator with compensated (Neumaier) summation.")
    setup = ["dPhiS_dx  = _ovec(T, 4)", "d2PhiS_dx2 = _ovec(T, 10)"]
    ret = """\
    return (PhiS_out,
            SVector{4,T}(dPhiS_dx[0], dPhiS_dx[1], dPhiS_dx[2], dPhiS_dx[3]),
            SVector{10,T}(d2PhiS_dx2[0], d2PhiS_dx2[1], d2PhiS_dx2[2], d2PhiS_dx2[3],
                          d2PhiS_dx2[4], d2PhiS_dx2[5], d2PhiS_dx2[6], d2PhiS_dx2[7],
                          d2PhiS_dx2[8], d2PhiS_dx2[9]),
            src_out)"""
    return build("calc", doc, "", "", setup, ret, "void effsource_ctx_calc(struct", comp=True)


def gen_phis_m():
    doc = ("m-mode singular field (complex), including phase exp(-i m phi_p).\n"
           "`comp=true` uses compensated summation for the numerator.")
    setup = ["PhiS = _ovec(T, 2)"]
    ret = "    return Complex(PhiS[0], PhiS[1])"
    return build("phis_m", doc, "m::Integer, ", "m, ", setup, ret,
                 "void effsource_ctx_PhiS_m(struct", comp=True)


def gen_calc_m():
    doc = ("m-mode singular field, derivatives and effective source (all complex).\n"
           "Returns `(PhiS, dPhiS_dx[4], d2PhiS_dx2[10], src)`.\n"
           "`comp=true` uses compensated summation for the numerator and source.")
    setup = ["PhiS       = _ovec(T, 2)", "dPhiS_dx   = _ovec(T, 8)",
             "d2PhiS_dx2 = _ovec(T, 20)", "src        = _ovec(T, 2)"]
    ret = """\
    return (Complex(PhiS[0], PhiS[1]),
            SVector{4,Complex{T}}(Complex(dPhiS_dx[0], dPhiS_dx[1]),
                                  Complex(dPhiS_dx[2], dPhiS_dx[3]),
                                  Complex(dPhiS_dx[4], dPhiS_dx[5]),
                                  Complex(dPhiS_dx[6], dPhiS_dx[7])),
            SVector{10,Complex{T}}(Complex(d2PhiS_dx2[0], d2PhiS_dx2[1]),
                                   Complex(d2PhiS_dx2[2], d2PhiS_dx2[3]),
                                   Complex(d2PhiS_dx2[4], d2PhiS_dx2[5]),
                                   Complex(d2PhiS_dx2[6], d2PhiS_dx2[7]),
                                   Complex(d2PhiS_dx2[8], d2PhiS_dx2[9]),
                                   Complex(d2PhiS_dx2[10], d2PhiS_dx2[11]),
                                   Complex(d2PhiS_dx2[12], d2PhiS_dx2[13]),
                                   Complex(d2PhiS_dx2[14], d2PhiS_dx2[15]),
                                   Complex(d2PhiS_dx2[16], d2PhiS_dx2[17]),
                                   Complex(d2PhiS_dx2[18], d2PhiS_dx2[19])),
            Complex(src[0], src[1]))"""
    return build("calc_m", doc, "m::Integer, ", "m, ", setup, ret,
                 "void effsource_ctx_calc_m(struct", comp=True)


header = """\
# ============================================================================
# Auto-generated by tools/generate_functions.py from kerr-circular-ctx.c.
# Do not edit by hand -- re-run the generator instead.
# ============================================================================
"""

OUT.write_text(
    header + "\n"
    + gen_phis() + "\n"
    + gen_calc() + "\n"
    + gen_phis_m() + "\n"
    + gen_calc_m() + "\n"
)
print(f"Wrote {OUT}")
