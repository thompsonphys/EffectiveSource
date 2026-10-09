#!/usr/bin/env python3
"""Self-tests for the transpiler constructs used by the equatorial kernels."""
from c2jl import convert_expr
from transpile import transpile_body


def jl(c, comp=True):
    return [line.strip() for line in transpile_body(c, comp=comp)]


def test_minus_assign():
    assert jl("x[0] -= y[0] * rt;") == ["x[0] -= y[0] * rt"]


def test_plus_assign_with_negative_rhs():
    assert jl("x[1] += -y[1] * rtt;") == ["x[1] += -y[1] * rtt"]


def test_zero_based_loop_unchanged():
    assert jl("for (int i = 0; i < 2; i++) s += v[i];", comp=False) == [
        "for i in 0:(2)-1", "s += v[i]", "end"]


def test_general_inclusive_loop():
    assert jl("for (int k = max(j - i, 0); k <= m + 2 + j; k++) s += C[k];", comp=False) == [
        "for k in max(j - i, 0):(m + 2 + j)", "s += C[k]", "end"]


def test_ellke_out_pointer_call():
    assert jl("effsource_ellKE(ctx, C1, &ellK, &ellE);") == ["(ellK, ellE) = _ellKE(ctx, C1)"]


def test_imei_table():
    assert convert_expr("ImEI[m][i][j][k] * fac2") == "imei(T,m,i,j,k) * fac2"


def test_comp_false_keeps_plain_source():
    assert jl("*src = -((a - b))/(c);", comp=False) == ["src_out = -((a - b))/(c)"]


def test_comp_true_compensates_source():
    assert "_csum(cm," in jl("*src = -((a - b))/(c);", comp=True)[0]


def test_comp_false_skips_accumulator():
    assert jl("num = 0; for (int i = 0; i < 2; i++) num += x[i];", comp=False) == [
        "num = 0", "for i in 0:(2)-1", "num += x[i]", "end"]


def test_comp_true_compensates_accumulator():
    assert any("_cadd(cm" in line for line in jl("num = 0; for (int i = 0; i < 2; i++) num += x[i];"))


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
    print("transpile self-test OK")
