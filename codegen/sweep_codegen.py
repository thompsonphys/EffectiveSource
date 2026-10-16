"""Generate the azimuthal-sweep form of effsource_equatorial_ctx_calc.

The sweep_calc function evaluates the polynomials of the calc body from coefficients in powers
of sin(dphib/2) and sin(dphib). The sweep_src function evaluates src alone from the exact form
s2^(11/2) src = sum_j (P_j + dQ cQ Q_j) dQ^(2j), whose coefficients the prepare function
computes once per field point and particle context.

Usage: python sweep_codegen.py CTX_SOURCE OUT_HEADER OUT_SOURCE
"""

import re
import sys
from math import comb

import sympy as sp
from sympy.printing.c import C99CodePrinter

CALC_NAME = "effsource_equatorial_ctx_calc"
SPLIT_MARKER = "/* Convert partial derivatives to total derivatives */"
SWEEP_STRUCT = "effsource_equatorial_sweep"
PREPARE_NAME = "effsource_equatorial_ctx_sweep_prepare"
PREPARE_DOUBLE_NAME = "effsource_equatorial_ctx_sweep_prepare_double"
PREPARE_DD_NAME = "effsource_equatorial_ctx_sweep_prepare_dd"
DD_THRESHOLD = "1e-2"
SWEEP_CALC_NAME = "effsource_equatorial_ctx_sweep_calc"
SWEEP_SRC_NAME = "effsource_equatorial_ctx_sweep_src"
THETA_TRIG = {"sin(theta)": "sinth", "sin(2 * theta)": "sin2th", "cos(2 * theta)": "cos2th"}
POWER_ALIAS = re.compile(r"^(dr|dtheta|dQ)(\d+)$")
IDENTIFIER = re.compile(r"[A-Za-z_]\w*")
CTX_FIELD = re.compile(r"^(A|dAdt|d2Adt2)\d+$")


def extract_function_body(text, name):
    """Return the body of the C function `name`, without the enclosing braces."""
    start = text.index(f"void {name}(")
    brace = text.index("{", start)
    end = text.index("\n}", brace)
    return text[brace + 1:end]


def strip_comments(text):
    """Remove C block comments."""
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def split_top_level(text, separator):
    """Split text at separators that are outside parentheses and brackets."""
    parts, depth, current = [], 0, []
    for char in text:
        if char in "([":
            depth += 1
        elif char in ")]":
            depth -= 1
        if char == separator and depth == 0:
            parts.append("".join(current))
            current = []
        else:
            current.append(char)
    parts.append("".join(current))
    return parts


def parse_statements(body):
    """Parse straight-line C statements into declarations and assignments."""
    statements = []
    for raw in split_top_level(strip_comments(body), ";"):
        text = " ".join(raw.split())
        if not text:
            continue
        declaration = re.match(r"^((?:const\s+)?(?:struct\s+\w+|double))\s+(.*)$", text)
        if declaration:
            declarators = []
            for item in split_top_level(declaration.group(2), ","):
                name, _, value = item.partition("=")
                declarators.append([name.strip(), value.strip() or None])
            statements.append({"kind": "decl", "type": declaration.group(1),
                               "declarators": declarators})
            continue
        assignment = re.match(r"^(\*?[A-Za-z_]\w*(?:\[\d+\])?)\s*=\s*(.*)$", text)
        if not assignment:
            raise ValueError(f"cannot parse statement: {text}")
        statements.append({"kind": "assign", "lhs": assignment.group(1),
                           "rhs": assignment.group(2)})
    return statements


def render(statement):
    """Render a parsed statement as one line of C."""
    if statement["kind"] == "assign":
        return f"  {statement['lhs']} = {statement['rhs']};"
    items = [name if value is None else f"{name} = {value}"
             for name, value in statement["declarators"]]
    return f"  {statement['type']} {', '.join(items)};"


def to_sympy(expression, symbols):
    """Parse a C arithmetic expression into sympy with exact rational literals."""
    names = set(IDENTIFIER.findall(expression))
    for name in names:
        symbols.setdefault(name, sp.Symbol(name))
    return sp.sympify(expression, locals=symbols, rational=True)


def power_aliases(statements):
    """Map dr2, dtheta4, dQ8 and similar names to powers of their base symbol."""
    aliases = {}
    for statement in statements:
        if statement["kind"] != "decl":
            continue
        for name, value in statement["declarators"]:
            match = POWER_ALIAS.match(name)
            if match and value is not None:
                aliases[sp.Symbol(name)] = sp.Symbol(match.group(1)) ** int(match.group(2))
    return aliases


def collect_polynomials(head_statements):
    """Collect each A-polynomial of the calc body by powers of dQ and dR."""
    aliases = power_aliases(head_statements)
    symbols = {}
    dQ, dR, dr, dtheta = (sp.Symbol(name) for name in ("dQ", "dR", "dr", "dtheta"))
    polynomials = {}
    for statement in head_statements:
        if statement["kind"] != "assign" or statement["lhs"].startswith("*"):
            continue
        expression = to_sympy(statement["rhs"], symbols).xreplace(aliases)
        if not expression.has(dQ) and not expression.has(dR):
            continue
        free = {symbol.name for symbol in expression.free_symbols}
        unexpected = {name for name in free
                      if name not in ("dQ", "dR", "dr", "dtheta") and not CTX_FIELD.match(name)}
        if unexpected:
            raise ValueError(f"{statement['lhs']} depends on {sorted(unexpected)}")
        polynomial = sp.Poly(sp.expand(expression), dQ, dR)
        if polynomial.degree(dR) > 1:
            raise ValueError(f"{statement['lhs']} is not linear in dR")
        parts = {}
        for (power_q, power_r), coefficient in polynomial.terms():
            parts.setdefault(power_r, {})[power_q] = coefficient
        polynomials[statement["lhs"]] = parts
    return polynomials, (dr, dtheta)


def horner_layout(terms):
    """Choose a Horner variable and return (variable, prefactor, coefficients by index)."""
    powers = sorted(terms)
    if all(power % 2 == 0 for power in powers):
        return "dQ2", None, {power // 2: terms[power] for power in powers}
    if all(power % 2 == 1 for power in powers):
        return "dQ2", "dQ", {(power - 1) // 2: terms[power] for power in powers}
    return "dQ", None, dict(terms)


def build_layouts(polynomials):
    """Return the struct fields and Horner layouts of every polynomial part."""
    layouts = []
    for name, parts in polynomials.items():
        for power_r in sorted(parts):
            variable, prefactor, coefficients = horner_layout(parts[power_r])
            size = max(coefficients) + 1
            field = f"{name}_R{power_r}"
            layouts.append({"name": name, "power_r": power_r, "field": field,
                            "variable": variable, "prefactor": prefactor,
                            "size": size, "coefficients": coefficients})
    return layouts


def horner_text(layout):
    """Return the C Horner expression of one polynomial part."""
    field = f"sweep->{layout['field']}"
    expression = f"{field}[{layout['size'] - 1}]"
    for index in range(layout["size"] - 2, -1, -1):
        expression = f"{field}[{index}] + {layout['variable']} * ({expression})"
    if layout["prefactor"]:
        expression = f"{layout['prefactor']} * ({expression})"
    return f"({expression})"


class DoublePrinter(C99CodePrinter):
    """C printer for double code: exact rational literals and integer powers as products."""

    def _print_Rational(self, expr):
        """Print p/q as a double division."""
        return f"{expr.p}.0/{expr.q}.0"

    def _print_Pow(self, expr):
        """Print integer powers as repeated products."""
        base, exponent = expr.as_base_exp()
        if exponent.is_Integer and 0 < abs(int(exponent)) <= 16:
            product = "*".join([f"({self._print(base)})"] * abs(int(exponent)))
            return f"({product})" if exponent > 0 else f"(1.0/({product}))"
        raise ValueError(f"unsupported power {expr}")


def dd_code(expression, names):
    """Return C code that evaluates expression in double-double arithmetic."""
    if expression.is_Symbol:
        return names.get(expression.name, expression.name)
    if expression.is_Integer:
        return f"dd_d({int(expression)}.0)"
    if expression.is_Rational:
        if expression.q & (expression.q - 1) == 0:
            return f"dd_d({expression.p}.0 / {expression.q}.0)"
        return f"dd_div(dd_d({expression.p}.0), dd_d({expression.q}.0))"
    if expression.is_Add:
        code = None
        for term in expression.args:
            coefficient, _ = term.as_coeff_Mul()
            negative = coefficient.is_negative
            piece = dd_code(-term if negative else term, names)
            if code is None:
                code = f"dd_neg({piece})" if negative else piece
            else:
                code = f"dd_{'sub' if negative else 'add'}({code}, {piece})"
        return code
    if expression.is_Mul:
        coefficient, rest = expression.as_coeff_Mul()
        code = None
        for factor in sp.Mul.make_args(rest):
            piece = dd_code(factor, names)
            code = piece if code is None else f"dd_mul({code}, {piece})"
        if coefficient == 1:
            return code
        if coefficient == -1:
            return f"dd_neg({code})"
        if coefficient.is_Integer or (coefficient.is_Rational and coefficient.q & (coefficient.q - 1) == 0):
            return f"dd_mul_d({code}, {coefficient.p}.0 / {coefficient.q}.0)"
        return f"dd_mul({dd_code(coefficient, names)}, {code})"
    if expression.is_Pow and expression.exp.is_Integer:
        exponent = int(expression.exp)
        base = dd_code(expression.base, names)
        code = base
        for _ in range(abs(exponent) - 1):
            code = f"dd_mul({code}, {base})"
        return code if exponent > 0 else f"dd_div(dd_d(1.0), {code})"
    raise ValueError(f"cannot write {expression} in double-double arithmetic")


def coefficient_symbol(layout, index):
    """Return the symbol of one coefficient of a polynomial part."""
    return sp.Symbol(f"{layout['field']}_{index}")


def evaluate(value, env, trig):
    """Parse a C expression with the symbolic environment and the azimuthal substitutions."""
    for pattern, field in THETA_TRIG.items():
        value = value.replace(pattern, field)
    local = {name: env.get(name, sp.Symbol(name)) for name in set(IDENTIFIER.findall(value))}
    local.update({"sin": sp.sin, "cos": sp.cos, "sqrt": sp.sqrt})
    return sp.sympify(value, locals=local, rational=True).xreplace(trig)


def symbolic_numerator(statements, polynomials):
    """Execute the calc body symbolically and return N = s2^(11/2) src and the symbols it uses.

    sqrt(s2) is replaced by s2^-5, so that every derivative of PhiS arrives multiplied by the
    power of s2 that clears its denominator, and src becomes N directly.
    """
    dphib, dQ, cQ, s2 = sp.symbols("dphib dQ cQ s2")
    trig = {sp.sin(dphib / 2): dQ, sp.cos(dphib / 2): cQ,
            sp.sin(dphib): 2 * dQ * cQ, sp.cos(dphib): 1 - 2 * dQ**2}
    env, s2_definition, numerator = {}, None, None
    for statement in statements:
        if statement["kind"] == "decl":
            pairs = [(name, value) for name, value in statement["declarators"] if value is not None]
        else:
            pairs = [(statement["lhs"], statement["rhs"])]
        for name, value in pairs:
            if name in polynomials and not statement.get("tail"):
                env[name] = sp.Symbol(name)
            elif name == "*src":
                numerator = evaluate(value, env, trig)
            elif name.startswith("*") or "[" in name:
                continue
            elif value.startswith("ctx->") or value.startswith("xp.") \
                    or name in ("r", "theta", "phi", "dr", "dtheta"):
                env[name] = sp.Symbol(name)
            elif name == "dphib":
                env[name] = dphib
            elif name == "s2":
                s2_definition = evaluate(value, env, trig)
                env[name] = s2
            elif name == "sqrts2":
                env[name] = s2**-5
            else:
                env[name] = evaluate(value, env, trig)
    return numerator, s2_definition, (dQ, cQ, s2)


def polynomial_form(expression, gens, families):
    """Return expression as {(power dQ, power cQ, power s2, family): scalar coefficient}.

    Scalar coefficients are kept unexpanded. The expression must be polynomial in dQ, cQ and
    s2 and linear in the family symbols.
    """
    dQ, cQ, s2 = gens
    special = {dQ: (1, 0, 0, None), cQ: (0, 1, 0, None), s2: (0, 0, 1, None)}
    special.update({family: (0, 0, 0, family.name) for family in families})
    tracked = set(special)
    cache = {}

    def add(left, right):
        result = dict(left)
        for key, value in right.items():
            result[key] = result[key] + value if key in result else value
        return result

    def multiply(left, right):
        result = {}
        for (q1, c1, s1, f1), v1 in left.items():
            for (q2, c2, s2_, f2), v2 in right.items():
                if f1 is not None and f2 is not None:
                    raise ValueError("s2^(11/2) src is not linear in the A polynomials")
                key = (q1 + q2, c1 + c2, s1 + s2_, f1 or f2)
                product = v1 * v2
                result[key] = result[key] + product if key in result else product
        return result

    def walk(node):
        if node in cache:
            return cache[node]
        if not (node.free_symbols & tracked):
            result = {(0, 0, 0, None): node}
        elif node in special:
            result = {special[node]: sp.Integer(1)}
        elif node.is_Add:
            result = {}
            for argument in node.args:
                result = add(result, walk(argument))
        elif node.is_Mul:
            result = {(0, 0, 0, None): sp.Integer(1)}
            for argument in node.args:
                result = multiply(result, walk(argument))
        elif node.is_Pow and node.exp.is_Integer and node.exp > 0:
            base = walk(node.base)
            result = {(0, 0, 0, None): sp.Integer(1)}
            for _ in range(int(node.exp)):
                result = multiply(result, base)
        else:
            raise ValueError(f"cannot write {node} as a polynomial in dQ, cQ and s2")
        cache[node] = result
        return result

    return walk(expression)


def rational_src_coefficients(statements, polynomials, layouts):
    """Return the weights and the sparse terms of P_j, Q_j in s2^(11/2) src = sum_j (P_j + dQ cQ Q_j) dQ^(2j).

    Each term is (integer factor, coefficient symbol, weight symbol). The weights are kept in the
    factored form of the calc body.
    """
    numerator, s2_definition, gens = symbolic_numerator(statements, polynomials)
    dQ, cQ, s2 = gens
    families = [sp.Symbol(name) for name in polynomials]
    form = polynomial_form(numerator, gens, families)
    alpha, beta = sp.symbols("alpha_s2 beta")
    alpha20, alpha02, dr, dtheta = sp.symbols("alpha20 alpha02 dr dtheta")
    if sp.expand(s2_definition - (alpha20 * dr**2 + alpha02 * dtheta**2 + beta * dQ**2)) != 0:
        raise ValueError(f"unexpected s2 = {s2_definition}")
    grouped = {}
    for (power_q, power_c, power_s, family), value in form.items():
        if family is None:
            if value != 0:
                raise ValueError("s2^(11/2) src has a term without the A polynomials")
            continue
        for m in range(power_s + 1):
            factor = comb(power_s, m) * alpha ** (power_s - m) * beta**m
            for n in range(power_c // 2 + 1):
                sign_factor = comb(power_c // 2, n) * (-1) ** n
                key = (family, power_q + 2 * m + 2 * n, power_c % 2)
                term = sign_factor * factor * value
                grouped[key] = grouped[key] + term if key in grouped else term
    weights, weight_symbols = [], {}
    for (family, power_q, power_c), value in grouped.items():
        symbol = sp.Symbol(f"w{len(weights)}")
        weights.append((symbol, value))
        weight_symbols[(family, (power_q, power_c))] = symbol
    family_terms = {name: [] for name in polynomials}
    for layout in layouts:
        for index in range(layout["size"]):
            power_q = 2 * index + (1 if layout["prefactor"] else 0)
            symbol = coefficient_symbol(layout, index)
            if layout["power_r"] == 1:
                family_terms[layout["name"]].append((power_q + 1, 1, 2, symbol))
            else:
                family_terms[layout["name"]].append((power_q, 0, 1, symbol))
    combined = {}
    for (name, (weight_q, weight_c)), weight in weight_symbols.items():
        for family_q, family_c, factor, symbol in family_terms[name]:
            power_q, power_c = weight_q + family_q, weight_c + family_c
            for m in range(power_c // 2 + 1):
                key = (power_q + 2 * m, power_c % 2)
                value = factor * comb(power_c // 2, m) * (-1) ** m
                bucket = combined.setdefault(key, {})
                bucket[(symbol, weight)] = bucket.get((symbol, weight), 0) + value
    for (power_q, power_c), bucket in combined.items():
        if any(bucket.values()) and (power_q + power_c) % 2:
            raise ValueError(f"s2^(11/2) src has the wrong parity at dQ^{power_q} cQ^{power_c}")
    size_p = max(power_q for power_q, power_c in combined if power_c == 0) // 2 + 1
    size_q = max(power_q for power_q, power_c in combined if power_c == 1) // 2 + 1

    def terms(key):
        return [(value, symbol, weight) for (symbol, weight), value in combined.get(key, {}).items()
                if value != 0]

    p_terms = [terms((2 * j, 0)) for j in range(size_p)]
    q_terms = [terms((2 * j + 1, 1)) for j in range(size_q)]
    return weights, p_terms, q_terms


def polynomial_text(name, layouts):
    """Return the C expression that evaluates the polynomial `name` from its parts."""
    pieces = []
    for layout in layouts:
        if layout["name"] != name:
            continue
        piece = horner_text(layout)
        pieces.append(piece if layout["power_r"] == 0 else f"dR * {piece}")
    return " + ".join(pieces) if pieces else "0.0"


def prepare_source(layouts, gens, ctx_scalars, weights, p_terms, q_terms, arithmetic):
    """Return the C body of a prepare function in "double" or "double-double" arithmetic."""
    dr, dtheta = gens
    names, expressions = [], []
    for layout in layouts:
        for index in range(layout["size"]):
            coefficient = layout["coefficients"].get(index, sp.Integer(0))
            if coefficient != 0:
                coefficient = sp.horner(sp.Poly(coefficient, dr, dtheta))
            names.append((layout["field"], index))
            expressions.append(coefficient)
    replacements, reduced = sp.cse(expressions, symbols=sp.numbered_symbols("t"))
    weight_replacements, weight_reduced = sp.cse([value for _, value in weights],
                                                 symbols=sp.numbered_symbols("v"))
    used = set()
    for expression in ([value for _, value in replacements] + reduced
                       + [value for _, value in weight_replacements] + weight_reduced):
        used |= {symbol.name for symbol in expression.free_symbols}
    scalars = sorted(name for name in used | {"alpha20", "alpha02"}
                     if CTX_FIELD.match(name) or name in ctx_scalars)
    lines = [
        "  sweep->r = r;",
        "  sweep->theta = theta;",
        "  sweep->sinth = sin(theta);",
        "  sweep->sin2th = sin(2 * theta);",
        "  sweep->cos2th = cos(2 * theta);",
    ]
    if arithmetic == "double":
        printer = DoublePrinter()
        kind = "double"

        def code(expression):
            return printer.doprint(expression)

        def to_double(value):
            return value

        def store_split(field, field_lo, value):
            return [f"  sweep->{field} = {value};", f"  sweep->{field_lo} = 0.0;"]

        lines += [
            "  const double dr = r - ctx->xp.r;",
            "  const double dtheta = theta - ctx->xp.theta;",
            "  const double sinth = sweep->sinth;",
            "  const double sin2th = sweep->sin2th;",
            "  const double cos2th = sweep->cos2th;",
        ]
        lines += [f"  const double {name} = ctx->{name};" for name in scalars]
    else:
        kind = "effsource_dd"
        rename = {"r": "r_dd"}

        def code(expression):
            return dd_code(expression, rename)

        def to_double(value):
            return f"dd_to_double({value})"

        def store_split(field, field_lo, value):
            return [f"  sweep->{field} = {value}.hi;", f"  sweep->{field_lo} = {value}.lo;"]

        lines += [
            "  const effsource_dd r_dd = dd_d(r);",
            "  const effsource_dd dr = dd_two_sum(r, -ctx->xp.r);",
            "  const effsource_dd dtheta = dd_two_sum(theta, -ctx->xp.theta);",
            "  const effsource_dd half_pi = {1.5707963267948966, 6.123233995736766e-17};",
            "  effsource_dd sin_delta, cos_delta;",
            "  dd_sincos(dd_sub(dd_d(theta), half_pi), &sin_delta, &cos_delta);",
            "  const effsource_dd sinth = cos_delta;",
            "  const effsource_dd sin2th = dd_neg(dd_mul_d(dd_mul(cos_delta, sin_delta), 2.0));",
            "  const effsource_dd cos2th = dd_sub(dd_d(1.0), dd_mul_d(dd_mul(cos_delta, cos_delta), 2.0));",
        ]
        lines += [f"  const effsource_dd {name} = dd_d(ctx->{name});" for name in scalars]
    for symbol, value in replacements:
        lines.append(f"  const {kind} {symbol} = {code(value)};")
    for (field, index), value in zip(names, reduced):
        lines.append(f"  const {kind} {field}_{index} = {code(value)};")
        lines.append(f"  sweep->{field}[{index}] = {to_double(f'{field}_{index}')};")
    alpha = sp.Symbol("alpha20") * dr**2 + sp.Symbol("alpha02") * dtheta**2
    lines.append(f"  const {kind} alpha_s2 = {code(alpha)};")
    for symbol, value in weight_replacements:
        lines.append(f"  const {kind} {symbol} = {code(value)};")
    for (symbol, _), value in zip(weights, weight_reduced):
        lines.append(f"  const {kind} {symbol} = {code(value)};")
    for field, all_terms in (("src_P", p_terms), ("src_Q", q_terms)):
        for index, terms in enumerate(all_terms):
            total = sp.Add(*[factor * symbol * weight for factor, symbol, weight in terms])
            lines.append(f"  const {kind} {field}_{index} = {code(total)};")
            lines += store_split(f"{field}[{index}]", f"{field}_lo[{index}]", f"{field}_{index}")
    lines += store_split("src_alpha", "src_alpha_lo", "alpha_s2")
    return "\n".join(lines)


def sweep_src_body(size_p, size_q):
    """Return the C body of the rational src evaluation."""
    def horner(field, size):
        expression = f"sweep->{field}[{size - 1}]"
        for index in range(size - 2, -1, -1):
            expression = f"sweep->{field}[{index}] + u * ({expression})"
        return expression
    return f"""  const double dr = sweep->r - ctx->xp.r;
  const double dphib = (phi_field - ctx->xp.phi) - ctx->c * dr;
  const double dQ = sin(0.5 * dphib);
  const double cQ = cos(0.5 * dphib);
  const double u = dQ * dQ;
  const double s2 = sweep->src_alpha + ctx->beta * u;
  const double s2_2 = s2 * s2;
  const double P = {horner("src_P", size_p)};
  const double Q = {horner("src_Q", size_q)};
  return (P + dQ * cQ * Q) / (s2_2 * s2_2 * s2 * sqrt(s2));"""


def referenced_names(statements):
    """Return every identifier read by the statements."""
    names = set()
    for statement in statements:
        if statement["kind"] == "assign":
            names |= set(IDENTIFIER.findall(statement["rhs"]))
            names |= set(IDENTIFIER.findall(statement["lhs"].lstrip("*").split("[")[0]))
        else:
            for _, value in statement["declarators"]:
                if value is not None:
                    names |= set(IDENTIFIER.findall(value))
    return names


def prune_unused(statements):
    """Remove initialised declarations whose names are never read."""
    while True:
        changed = False
        for index, statement in enumerate(statements):
            if statement["kind"] != "decl":
                continue
            others = statements[:index] + statements[index + 1:]
            reads = referenced_names(others)
            for declarator in list(statement["declarators"]):
                name, value = declarator
                own = set()
                for other_name, other_value in statement["declarators"]:
                    if other_name != name and other_value is not None:
                        own |= set(IDENTIFIER.findall(other_value))
                if value is not None and name not in reads and name not in own:
                    statement["declarators"].remove(declarator)
                    changed = True
        statements[:] = [statement for statement in statements
                         if statement["kind"] != "decl" or statement["declarators"]]
        if not changed:
            return statements


def sweep_calc_body(statements, polynomials, layouts):
    """Return the C body of the sweep evaluation function."""
    rewritten = []
    for statement in statements:
        statement = {key: (list(map(list, value)) if key == "declarators" else value)
                     for key, value in statement.items()}
        if statement["kind"] == "assign" and statement["lhs"] in polynomials \
                and not statement.get("tail"):
            statement["rhs"] = polynomial_text(statement["lhs"], layouts)
        if statement["kind"] == "decl":
            for declarator in statement["declarators"]:
                name, value = declarator
                if value == "x->r":
                    declarator[1] = "sweep->r"
                elif value == "x->theta":
                    declarator[1] = "sweep->theta"
                elif value == "x->phi":
                    declarator[1] = "phi_field"
                elif value in THETA_TRIG:
                    declarator[1] = f"sweep->{THETA_TRIG[value]}"
        rewritten.append(statement)
    for statement in rewritten:
        if statement["kind"] == "assign":
            for pattern, field in THETA_TRIG.items():
                statement["rhs"] = statement["rhs"].replace(pattern, f"sweep->{field}")
    rewritten = prune_unused(rewritten)
    return "\n".join(render(statement) for statement in rewritten)


def header_text(layouts, size_p, size_q):
    """Return the generated header."""
    fields = "\n".join(f"  double {layout['field']}[{layout['size']}];" for layout in layouts)
    fields += f"\n  double src_P[{size_p}], src_Q[{size_q}], src_alpha;"
    fields += f"\n  double src_P_lo[{size_p}], src_Q_lo[{size_q}], src_alpha_lo;"
    return f"""/* Generated by codegen/sweep_codegen.py from kerr-equatorial-ctx.c. Do not edit. */

#ifndef EFFSOURCE_EQUATORIAL_SWEEP_H
#define EFFSOURCE_EQUATORIAL_SWEEP_H

#include "effsource_equatorial.h"

/* Coefficients of {CALC_NAME} at one field point (r, theta) and one particle
   context, collected by powers of dQ = sin(dphib/2) and dR = sin(dphib). */
struct {SWEEP_STRUCT} {{
  double r, theta;
  double sinth, sin2th, cos2th;
{fields}
}};

/* Separation below which {PREPARE_NAME} uses double-double arithmetic, in units of
   C1 = (alpha20 dr^2 + alpha02 dtheta^2) / beta. */
#define EFFSOURCE_SWEEP_DD_C1 {DD_THRESHOLD}

/* Compute the coefficients for the field point (r, theta): with {PREPARE_DD_NAME}
   if C1 < EFFSOURCE_SWEEP_DD_C1 and with {PREPARE_DOUBLE_NAME} otherwise. */
void {PREPARE_NAME}(const struct effsource_equatorial_ctx * ctx, double r, double theta,
  struct {SWEEP_STRUCT} * sweep);

/* Compute the coefficients in double-double arithmetic and round them once to
   double. The src coefficients src_P, src_Q and src_alpha also keep their
   low parts in src_P_lo, src_Q_lo and src_alpha_lo, which only the analytic
   m-modes use. The src coefficients cancel strongly near the particle and a
   coefficient error repeats identically at every azimuth, so double precision
   loses digits there. With this function the src of {SWEEP_SRC_NAME} is
   correct to about 1e-15 relative to its largest value over the azimuths. */
void {PREPARE_DD_NAME}(const struct effsource_equatorial_ctx * ctx, double r, double theta,
  struct {SWEEP_STRUCT} * sweep);

/* Compute the coefficients in double arithmetic: about 20 times faster than
   {PREPARE_DD_NAME}, and as accurate as {CALC_NAME} within a factor of a few.
   The low parts src_P_lo, src_Q_lo and src_alpha_lo are zero. */
void {PREPARE_DOUBLE_NAME}(const struct effsource_equatorial_ctx * ctx, double r, double theta,
  struct {SWEEP_STRUCT} * sweep);

/* Same outputs as {CALC_NAME} at the field point (sweep->r, sweep->theta, phi_field). */
void {SWEEP_CALC_NAME}(const struct effsource_equatorial_ctx * ctx,
  const struct {SWEEP_STRUCT} * sweep, double phi_field,
  double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);

/* The src output of {CALC_NAME} at (sweep->r, sweep->theta, phi_field), from
   s2^(11/2) src = sum_j (src_P[j] + dQ cQ src_Q[j]) dQ^(2j), with
   dQ = sin(dphib/2), cQ = cos(dphib/2) and s2 = src_alpha + beta dQ^2. */
double {SWEEP_SRC_NAME}(const struct effsource_equatorial_ctx * ctx,
  const struct {SWEEP_STRUCT} * sweep, double phi_field);

#endif /* EFFSOURCE_EQUATORIAL_SWEEP_H */
"""


def source_text(header_name, prepare_body, prepare_double_body, calc_body, src_body):
    """Return the generated C source."""
    return f"""/* Generated by codegen/sweep_codegen.py from kerr-equatorial-ctx.c. Do not edit. */

#include <math.h>
#include "effsource_dd.h"
#include "{header_name}"

void {PREPARE_DD_NAME}(const struct effsource_equatorial_ctx *ctx, double r, double theta,
                       struct {SWEEP_STRUCT} *sweep)
{{
{prepare_body}
}}

void {PREPARE_DOUBLE_NAME}(const struct effsource_equatorial_ctx *ctx, double r, double theta,
                           struct {SWEEP_STRUCT} *sweep)
{{
{prepare_double_body}
}}

void {PREPARE_NAME}(const struct effsource_equatorial_ctx *ctx, double r, double theta,
                    struct {SWEEP_STRUCT} *sweep)
{{
  const double dr = r - ctx->xp.r;
  const double dtheta = theta - ctx->xp.theta;
  const double C1 = (ctx->alpha20 * dr * dr + ctx->alpha02 * dtheta * dtheta) / ctx->beta;
  if (C1 < EFFSOURCE_SWEEP_DD_C1)
    {PREPARE_DD_NAME}(ctx, r, theta, sweep);
  else
    {PREPARE_DOUBLE_NAME}(ctx, r, theta, sweep);
}}

void {SWEEP_CALC_NAME}(const struct effsource_equatorial_ctx *ctx,
                       const struct {SWEEP_STRUCT} *sweep, double phi_field,
                       double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src)
{{
{calc_body}
}}

double {SWEEP_SRC_NAME}(const struct effsource_equatorial_ctx *ctx,
                        const struct {SWEEP_STRUCT} *sweep, double phi_field)
{{
{src_body}
}}
"""


def main():
    """Read the context source and write the sweep header and source."""
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    source_path, header_path, out_path = sys.argv[1:]
    text = open(source_path).read()
    body = extract_function_body(text, CALC_NAME)
    head_raw, _, tail_raw = body.partition(SPLIT_MARKER)
    if not tail_raw:
        raise ValueError("split marker not found in the calc body")
    head = parse_statements(head_raw)
    tail = parse_statements(tail_raw)
    for statement in tail:
        statement["tail"] = True
    polynomials, gens = collect_polynomials(head)
    layouts = build_layouts(polynomials)
    ctx_scalars = {name for statement in head if statement["kind"] == "decl"
                   for name, value in statement["declarators"]
                   if value == f"ctx->{name}" and not CTX_FIELD.match(name)}
    weights, p_terms, q_terms = rational_src_coefficients(head + tail, polynomials, layouts)
    prepare_body = prepare_source(layouts, gens, ctx_scalars, weights, p_terms, q_terms,
                                  "double-double")
    prepare_double_body = prepare_source(layouts, gens, ctx_scalars, weights, p_terms, q_terms,
                                         "double")
    calc_body = sweep_calc_body(head + tail, polynomials, layouts)
    src_body = sweep_src_body(len(p_terms), len(q_terms))
    header_name = header_path.rsplit("/", 1)[-1]
    with open(header_path, "w") as stream:
        stream.write(header_text(layouts, len(p_terms), len(q_terms)))
    with open(out_path, "w") as stream:
        stream.write(source_text(header_name, prepare_body, prepare_double_body, calc_body, src_body))
    for layout in layouts:
        print(f"{layout['field']:>14}: {layout['variable']:>3} x{layout['size']:<2}"
              f"{' * ' + layout['prefactor'] if layout['prefactor'] else ''}")
    print(f"rational src: P x{len(p_terms)}, Q x{len(q_terms)}, {len(weights)} weights, "
          f"{sum(len(terms) for terms in p_terms + q_terms)} products")


if __name__ == "__main__":
    main()
