#!/usr/bin/env python3
"""Build Coordinate and Offset Julia methods from a C kernel body."""
from transpile import transpile_body

RET_PHIS = "    return PhiS_out"
SETUP_CALC = ["dPhiS_dx  = _ovec(T, 4)", "d2PhiS_dx2 = _ovec(T, 10)"]
RET_CALC = """\
    return (PhiS_out,
            SVector{4,T}(dPhiS_dx[0], dPhiS_dx[1], dPhiS_dx[2], dPhiS_dx[3]),
            SVector{10,T}(d2PhiS_dx2[0], d2PhiS_dx2[1], d2PhiS_dx2[2], d2PhiS_dx2[3],
                          d2PhiS_dx2[4], d2PhiS_dx2[5], d2PhiS_dx2[6], d2PhiS_dx2[7],
                          d2PhiS_dx2[8], d2PhiS_dx2[9]),
            src_out)"""
SETUP_PHIS_M = ["PhiS = _ovec(T, 2)"]
RET_PHIS_M = "    return Complex(PhiS[0], PhiS[1])"
SETUP_CALC_M = ["PhiS       = _ovec(T, 2)", "dPhiS_dx   = _ovec(T, 8)",
                "d2PhiS_dx2 = _ovec(T, 20)", "src        = _ovec(T, 2)"]
RET_CALC_M = """\
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

_COORD_LHS = {"r", "theta", "phi", "rp", "thetap", "phip", "dr", "dtheta", "dphi"}
_COORD_RHS = {"x.r", "x.theta", "x.phi", "xp.r", "xp.theta", "xp.phi",
              "r - rp", "theta - thetap", "phi - phip"}


def _match_brace(s, open_idx):
    depth = 0
    for i in range(open_idx, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced braces")


def body_after(src, sig_marker):
    """Extract the { ... } body following the signature containing sig_marker."""
    start = src.index(sig_marker)
    brace = src.index("{", start)
    return src[brace + 1 : _match_brace(src, brace)]


def indent(lines, pad="    "):
    return "\n".join(pad + l if l else l for l in lines)


def _is_coord_line(line):
    s = line.strip()
    if "=" not in s:
        return False
    lhs, rhs = s.split("=", 1)
    return lhs.strip() in _COORD_LHS and rhs.strip() in _COORD_RHS


def offset_body(body):
    """Return (core, has_dphi): body with the coordinate preamble replaced by the offset preamble."""
    idxs = [i for i, l in enumerate(body) if _is_coord_line(l)]
    first = idxs[0]
    idxset = set(idxs)
    has_dphi = any(body[i].strip().startswith("dphi ") for i in idxs)
    core = [l for i, l in enumerate(body) if i not in idxset]
    pre = ["rp = xp.r", "thetap = xp.theta", "dr = off.dr", "dtheta = off.dtheta"]
    if has_dphi:
        pre.append("dphi = off.dphi")
    pre += ["r = rp + dr", "theta = thetap + dtheta"]
    return core[:first] + pre + core[first:], has_dphi


def build(src, ctx_type, name, doc, extra_args, extra_call, setup, retblock, sig_marker,
          comp, post=None):
    """Emit the Coordinate and Offset methods of one kernel sharing an identical core."""
    body = transpile_body(body_after(src, sig_marker), comp=comp)
    if post is not None:
        body = post(body)
    off, _ = offset_body(body)
    setup_block = ("\n" + indent(setup)) if setup else ""
    ea = extra_args
    eac = extra_call

    def core(argname, argtype, b):
        fname = f"_{name}" if comp else name
        cmarg = ", cm::Val" if comp else ""
        return (
            f"function {fname}(ctx::{ctx_type}{{T}}, {ea}{argname}::{argtype}{{T}}{cmarg}) where {{T}}"
            + setup_block + "\n" + indent(b) + "\n" + retblock + "\nend\n"
        )

    def wrapper(argname, argtype):
        return (
            f"{name}(ctx::{ctx_type}{{T}}, {ea}{argname}::{argtype}{{T}}; comp::Bool=false) where {{T}} = "
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
