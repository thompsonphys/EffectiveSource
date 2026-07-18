#!/usr/bin/env python3
"""C function-body -> Julia transpiler for the four kerr-circular-ctx.c kernels.

Handles exactly the constructs those functions use: local scalar/array decls,
initializer lists, element assignments, compound `+=`, braceless nested `for`
loops, an `if(m>20){...return}` guard, GSL elliptic-integral calls, and pointer
dereference outputs. Local arrays become 0-based OffsetArrays so C indices map
verbatim. Numeric genericity (type T) is handled by c2jl.convert_expr.
"""
import re
from c2jl import convert_expr, _match_paren, _split_top_comma, split_add_terms

# ---------------------------------------------------------------------------
# Preprocessing
# ---------------------------------------------------------------------------

def strip_comments(code):
    code = re.sub(r"/\*.*?\*/", " ", code, flags=re.DOTALL)
    code = re.sub(r"//[^\n]*", " ", code)
    return code


def convert_gsl(expr):
    """gsl_sf_ellint_{E,K}comp(arg, GSL_PREC_DOUBLE) -> _ellE/_ellK(arg)."""
    for cname, jname in (("gsl_sf_ellint_Ecomp", "_ellE"),
                         ("gsl_sf_ellint_Kcomp", "_ellK")):
        while True:
            idx = expr.find(cname + "(")
            if idx == -1:
                break
            open_idx = idx + len(cname)
            close_idx = _match_paren(expr, open_idx)
            args = _split_top_comma(expr[open_idx + 1 : close_idx])
            arg0 = args[0].strip()
            expr = expr[:idx] + f"{jname}({arg0})" + expr[close_idx + 1 :]
    return expr


def convert_rhs(expr):
    """Full expression conversion incl. GSL calls."""
    return convert_expr(convert_gsl(expr))


def compensate_src(rhs_c):
    """Rewrite the effective-source RHS `-((SUM))/((DEN))` so its numerator SUM
    is summed with the opt-in compensated helper: `-(_csum(cm, t1, t2, ...))/DEN`.
    Falls back to a plain conversion if the shape is unexpected."""
    s = rhs_c.strip()
    op = s.find("(")
    if not s.startswith("-") or op == -1:
        return convert_rhs(rhs_c)
    close = _match_paren(s, op)
    inner = s[op + 1:close].strip()           # drop the -( ... ) grouping layer
    if inner.startswith("(") and _match_paren(inner, 0) == len(inner) - 1:
        inner = inner[1:-1].strip()           # drop the sum's own paren layer
    terms = []
    for t in split_add_terms(inner):
        t = t.lstrip()
        if t.startswith("+"):
            t = t[1:].strip()                 # unary plus is implicit
        terms.append(convert_rhs(t))
    return f"{s[:op]}(_csum(cm, {', '.join(terms)})){s[close + 1:]}"


# ---------------------------------------------------------------------------
# Recursive-descent parser for the C body
# ---------------------------------------------------------------------------

class Parser:
    def __init__(self, s):
        self.s = s
        self.i = 0
        self.n = len(s)

    def _ws(self):
        while self.i < self.n and self.s[self.i] in " \t\r\n":
            self.i += 1

    def _kw(self, kw):
        j = self.i
        if self.s[j : j + len(kw)] == kw:
            nxt = self.s[j + len(kw) : j + len(kw) + 1]
            if nxt == "" or not (nxt.isalnum() or nxt == "_"):
                return True
        return False

    def parse_items(self):
        items = []
        while True:
            self._ws()
            if self.i >= self.n or self.s[self.i] == "}":
                return items
            items.append(self.parse_item())

    def parse_item(self):
        self._ws()
        if self._kw("for"):
            return self.parse_for()
        if self._kw("if"):
            return self.parse_if()
        if self.s[self.i] == "{":
            return self.parse_block()
        return self.parse_simple()

    def _paren_group(self):
        self._ws()
        assert self.s[self.i] == "(", f"expected ( at {self.i}: {self.s[self.i:self.i+20]!r}"
        close = _match_paren(self.s, self.i)
        inner = self.s[self.i + 1 : close]
        self.i = close + 1
        return inner

    def parse_for(self):
        self.i += 3
        header = self._paren_group()
        body = self.parse_item()
        return ("for", header, body)

    def parse_if(self):
        self.i += 2
        cond = self._paren_group()
        body = self.parse_item()
        return ("if", cond, body)

    def parse_block(self):
        self.i += 1  # consume '{'
        items = self.parse_items()
        self._ws()
        assert self.s[self.i] == "}", "expected }"
        self.i += 1
        return ("block", items)

    def parse_simple(self):
        depth = 0
        start = self.i
        while self.i < self.n:
            ch = self.s[self.i]
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == ";" and depth == 0:
                text = self.s[start : self.i].strip()
                self.i += 1
                return ("stmt", text)
            self.i += 1
        raise ValueError(f"unterminated statement: {self.s[start:start+40]!r}")


# ---------------------------------------------------------------------------
# Emission
# ---------------------------------------------------------------------------

_TYPE_RE = re.compile(r"^(const\s+)?(struct\s+coordinate|double|int)\b\s*")
_ARR_DECL_RE = re.compile(r"^(\w+)\[(\d+)\]$")
_ARR_INIT_RE = re.compile(r"^(\w+)\[(\d+)\]\s*=\s*\{(.*)\}$", re.DOTALL)


def emit_stmt(text, accum=None):
    """Return a list of Julia lines for one C simple-statement. If `accum` is the
    name of a compensated accumulator, a `NAME += EXPR` for it is emitted as a
    Neumaier step `NAME, NAME_c = _cadd(cm, NAME, NAME_c, EXPR)`."""
    text = text.replace("->", ".").strip()
    if text == "" or text == "return":
        # bare return only appears inside the m>20 guard -> handled in emit_if
        return ["#= return =#"]
    if text.startswith("printf"):
        return []  # diagnostic; dropped (guard emits an error instead)

    m = _TYPE_RE.match(text)
    if m:
        rest = text[m.end():]
        lines = []
        for item in _split_top_comma(rest):
            item = item.strip()
            mi = _ARR_INIT_RE.match(item)
            if mi:
                name, n, lst = mi.group(1), int(mi.group(2)), mi.group(3)
                elems = [convert_rhs(e) for e in _split_top_comma(lst)]
                lines.append(f"{name} = _ovec0(T[{', '.join(elems)}])")
                continue
            ma = _ARR_DECL_RE.match(item)
            if ma:
                lines.append(f"{ma.group(1)} = _ovec(T, {ma.group(2)})")
                continue
            if "=" in item:
                lhs, rhs = item.split("=", 1)
                lhs = lhs.strip()
                rhs_jl = compensate_src(rhs) if lhs == "effsrc" else convert_rhs(rhs)
                lines.append(f"{lhs} = {rhs_jl}")
            # bare scalar decl -> dropped
        return lines

    # not a declaration: assignment / compound assignment / pointer output
    if "+=" in text:
        lhs, rhs = text.split("+=", 1)
        lhs = lhs.strip()
        if accum is not None and lhs == accum:
            return [f"{lhs}, {lhs}_c = _cadd(cm, {lhs}, {lhs}_c, {convert_rhs(rhs)})"]
        return [f"{lhs} += {convert_rhs(rhs)}"]
    lhs, rhs = text.split("=", 1)
    lhs = lhs.strip()
    is_src = False
    if lhs.startswith("*"):
        base = lhs[1:].strip()
        lhs = base + "_out"              # *src -> src_out, *PhiS -> PhiS_out
        is_src = base == "src"
    rhs_jl = compensate_src(rhs) if is_src else convert_rhs(rhs)
    return [f"{lhs} = {rhs_jl}"]


def _for_header(header):
    """`int i=0; i<N; i++` -> ('i', 'N')."""
    parts = [p.strip() for p in header.split(";")]
    var = re.match(r"int\s+(\w+)\s*=", parts[0]).group(1)
    bound = re.match(r"\w+\s*<\s*(\w+)", parts[1]).group(1)
    return var, bound


def emit(node, indent, out, accum=None):
    pad = "    " * indent
    kind = node[0]
    if kind == "stmt":
        for line in emit_stmt(node[1], accum):
            out.append(pad + line)
    elif kind == "block":
        for it in node[1]:
            emit(it, indent, out, accum)
    elif kind == "for":
        var, bound = _for_header(node[1])
        out.append(f"{pad}for {var} in 0:({bound})-1")
        emit(node[2], indent + 1, out, accum)
        out.append(f"{pad}end")
    elif kind == "if":
        cond = node[1].replace("->", ".").strip()
        # the sole `if` in these kernels is the m>20 unsupported-mode guard
        out.append(f"{pad}if {convert_rhs(cond)}")
        out.append(f"{pad}    error(\"effsource: requested m-mode is not supported (m>20)\")")
        out.append(f"{pad}end")


# ---------------------------------------------------------------------------
# C1: opt-in compensated accumulation. Detect the exact C shape
#   NAME = 0;  for(...) for(...) for(...) NAME += EXPR;
# and emit Neumaier-compensated code (init + companion counter, _cadd step,
# finalize). The compensation itself is a runtime no-op unless the caller passes
# comp=true and the type is an IEEE float (see _cadd in EffSourceCircular.jl).
# ---------------------------------------------------------------------------

def _assign_parts(text):
    """(lhs, op, rhs) for a simple assignment, after stripping a leading decl
    keyword and pointer arrows. op is '=' or '+=' (else all None)."""
    t = text.replace("->", ".").strip()
    m = _TYPE_RE.match(t)
    if m:
        t = t[m.end():].strip()
    if "+=" in t:
        lhs, rhs = t.split("+=", 1)
        return lhs.strip(), "+=", rhs.strip()
    if "=" in t:
        lhs, rhs = t.split("=", 1)
        return lhs.strip(), "=", rhs.strip()
    return None, None, None


def _innermost_stmt(node):
    if node[0] == "for":
        return _innermost_stmt(node[2])
    if node[0] == "block":
        return _innermost_stmt(node[1][-1]) if node[1] else None
    return node if node[0] == "stmt" else None


def _accum_name(init_node, for_node):
    """If `init_node` is `NAME = 0` and `for_node` is a loop accumulating into
    NAME via `NAME += ...`, return NAME; else None."""
    if init_node[0] != "stmt" or for_node is None or for_node[0] != "for":
        return None
    lhs, op, rhs = _assign_parts(init_node[1])
    if op != "=" or rhs != "0":
        return None
    inner = _innermost_stmt(for_node)
    if inner is None:
        return None
    ilhs, iop, _ = _assign_parts(inner[1])
    return lhs if (iop == "+=" and ilhs == lhs) else None


def transpile_body(c_body):
    p = Parser(strip_comments(c_body))
    items = p.parse_items()
    out = []
    i, n = 0, len(items)
    while i < n:
        nxt = items[i + 1] if i + 1 < n else None
        name = _accum_name(items[i], nxt) if nxt is not None else None
        if name:
            out.append(f"    {name} = zero(T)")
            out.append(f"    {name}_c = zero(T)")
            emit(nxt, 1, out, accum=name)
            out.append(f"    {name} = {name} + {name}_c")
            i += 2
        else:
            emit(items[i], 1, out)
            i += 1
    return out
