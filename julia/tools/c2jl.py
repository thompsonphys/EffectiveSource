#!/usr/bin/env python3
"""Small, focused C->Julia expression/body transpiler for kerr-circular-ctx.c.

It only supports the constructs that actually occur in the ported functions:
  * pow(a, b)              -> (a)^(b)   (nested, paren-matched)
  * /* ... */ and // ...   -> # comments
  * `double`/`const double` local declarations -> plain assignments or dropped
  * `double name[N];`      -> OffsetArray of length N, 0-based (0:N-1)
  * `for(int i=0;i<N;i++)` -> `for i in 0:N-1`
  * `x->field`, `p.field`  -> `x.field`
  * `*ptr = expr;`         -> handled by caller (returned raw for output packing)

Local temporary arrays use OffsetArrays so C indices (A[0], C[k], ellip[i]) map
verbatim -- no error-prone index shifting. `ReEI[m][i][j][k]` -> `reei(m,i,j,k)`.

This module exposes helpers used by generate_*.py scripts.
"""
import re


def _match_paren(s, open_idx):
    """Given index of '(' return index of the matching ')'."""
    depth = 0
    for i in range(open_idx, len(s)):
        if s[i] == "(":
            depth += 1
        elif s[i] == ")":
            depth -= 1
            if depth == 0:
                return i
    raise ValueError("unbalanced parens")


def _split_top_comma(s):
    """Split on commas that are at paren/bracket depth 0."""
    parts, depth, last = [], 0, 0
    for i, ch in enumerate(s):
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif ch == "," and depth == 0:
            parts.append(s[last:i])
            last = i + 1
    parts.append(s[last:])
    return parts


def _pow_repl(a, b):
    """Build a Julia power expression for base `a`, exponent `b`.

    Integer exponents stay as `(a)^(n)` (type-preserving). Half-integer
    exponents (n.5) become `(a)^(n) * sqrt(a)` so the result stays in type T
    instead of promoting to Float64 via a float exponent. The whole thing is
    parenthesized so it composes correctly inside larger expressions.
    """
    bs = b.strip()
    # half-integer literal exponent, e.g. 3.5 or 2.5
    m = re.fullmatch(r"(\d+)\.5", bs)
    if m:
        n = int(m.group(1))
        return f"(({a})^({n}) * sqrt({a}))"
    return f"(({a})^({b}))"


def split_add_terms(expr):
    """Split `expr` at top-level (depth-0) binary + / - into signed terms, e.g.
    `a - b*c + d` -> ['a', '- b*c', '+ d']. Operators inside parens/brackets are
    not split. Each returned term keeps its leading sign (the first has none)."""
    terms, depth, start, i, n = [], 0, 0, 0, len(expr)
    while i < n:
        ch = expr[i]
        if ch in "([{":
            depth += 1
        elif ch in ")]}":
            depth -= 1
        elif depth == 0 and ch in "+-" and i > 0:
            j = i - 1
            while j >= 0 and expr[j] == " ":
                j -= 1
            prev = expr[j] if j >= 0 else ""
            # binary (subtraction/addition) only if preceded by an operand end;
            # otherwise it is a unary sign belonging to the following token
            if prev and (prev.isalnum() or prev in ")]_."):
                terms.append(expr[start:i].strip())
                start = i
        i += 1
    terms.append(expr[start:].strip())
    return terms


def convert_pow(expr):
    """Replace pow(a,b) with a Julia power, innermost via rightmost-first scan."""
    while True:
        idx = expr.rfind("pow(")
        if idx == -1:
            return expr
        open_idx = idx + 3  # position of '('
        close_idx = _match_paren(expr, open_idx)
        inner = expr[open_idx + 1 : close_idx]
        args = _split_top_comma(inner)
        if len(args) != 2:
            raise ValueError(f"pow with {len(args)} args: {inner!r}")
        a, b = args[0].strip(), args[1].strip()
        repl = _pow_repl(a, b)
        expr = expr[:idx] + repl + expr[close_idx + 1 :]


# Float literal not preceded by an identifier char or a dot (avoids matching
# the "008" in A008 or a decimal that is part of something else).
_FLOAT_RE = re.compile(r"(?<![A-Za-z0-9_.])(\d+\.\d*|\.\d+)([eE][+-]?\d+)?")


def wrap_float_literals(expr):
    """Wrap Float64 literals as T(...) so they don't promote a low-precision T
    up to Float64. Integer literals are left alone (they promote cleanly)."""
    def repl(m):
        lit = m.group(0)
        if lit.endswith("."):
            lit += "0"          # 3. -> 3.0
        elif lit.startswith("."):
            lit = "0" + lit      # .5 -> 0.5
        return f"T({lit})"
    return _FLOAT_RE.sub(repl, expr)


def convert_reei(expr):
    """ReEI[m][i][j][k] -> reei(T,m,i,j,k) (T so constants stay in-type)."""
    pat = re.compile(r"ReEI\[([^\]]*)\]\[([^\]]*)\]\[([^\]]*)\]\[([^\]]*)\]")
    return pat.sub(lambda m: f"reei(T,{m.group(1)},{m.group(2)},{m.group(3)},{m.group(4)})", expr)


def convert_expr(expr):
    """Convert a single C expression (no statements) to Julia."""
    expr = convert_reei(expr)
    expr = convert_pow(expr)
    expr = wrap_float_literals(expr)
    expr = re.sub(r"\bNAN\b", "T(NaN)", expr)   # C NAN -> typed NaN
    return expr.strip()


if __name__ == "__main__":
    # tiny self-test
    tests = [
        ("pow(a,2)", "((a)^(2))"),
        ("pow(pow(a,2),3)", "((((a)^(2)))^(3))"),
        ("64*pow(pow(L,2) + pow(r,2),3)", "64*((((L)^(2)) + ((r)^(2)))^(3))"),
        ("pow(r,3)+pow(a,2)", "((r)^(3))+((a)^(2))"),
    ]
    for src, want in tests:
        got = convert_pow(src)
        assert got == want, f"{src!r}: got {got!r} want {want!r}"
    assert convert_reei("ReEI[m][i][j][k]*x") == "reei(T,m,i,j,k)*x"
    # half-integer power decomposition
    assert convert_pow("pow(rho2,3.5)") == "((rho2)^(3) * sqrt(rho2))"
    assert convert_pow("A/pow(rho2,3.5)") == "A/((rho2)^(3) * sqrt(rho2))"
    assert convert_pow("pow(alpha+beta,2.5)") == "((alpha+beta)^(2) * sqrt(alpha+beta))"
    # float-literal wrapping, integers untouched
    assert wrap_float_literals("3.*x + 64*y") == "T(3.0)*x + 64*y"
    assert wrap_float_literals("0.5*d") == "T(0.5)*d"
    assert wrap_float_literals("A008 + 4.0*z") == "A008 + T(4.0)*z"
    assert wrap_float_literals("(x)^(8)") == "(x)^(8)"
    # term splitting: top-level +/- only, parens protected
    assert split_add_terms("2*a - b*c - d") == ["2*a", "- b*c", "- d"]
    assert split_add_terms("a*(x + y) - b") == ["a*(x + y)", "- b"]
    assert split_add_terms("f + g*(-2 + r)*h - k") == ["f", "+ g*(-2 + r)*h", "- k"]
    print("c2jl self-test OK")
