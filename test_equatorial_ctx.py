"""Test that the equatorial context-based API produces identical results to the global API."""
import math
import effsource_equatorial as es

# Test parameters
M = 1.0
a = 0.5
r_p = 10.0
theta_p = math.pi / 2.0
phi_p = math.pi / 3.0
t_p = 1.0

# Orbital parameters for an eccentric equatorial orbit
E = 0.95
L = 3.5
ur_p = 0.1

es.disable_gsl_error_handler()

# Set up particle
x_p = es.make_coordinate(t=t_p, r=r_p, theta=theta_p, phi=phi_p)

# --- Global API ---
es.effsource_init(M, a)
es.effsource_set_particle(x_p, E, L, ur_p)

# Field point slightly offset from particle
x = es.make_coordinate(t=0.0, r=r_p + 0.1, theta=theta_p + 0.01, phi=phi_p + 0.02)

global_PhiS = es.calc_PhiS(x)
global_calc = es.calc(x)
global_PhiS_m = es.calc_PhiS_m(2, x)
global_calc_m = es.calc_m(2, x)

# --- Context API ---
ctx = es.EffsourceEquatorialContext(M, a)
ctx.set_particle(x_p, E, L, ur_p)

ctx_PhiS = ctx.calc_PhiS(x)
ctx_calc = ctx.calc(x)
ctx_PhiS_m = ctx.calc_PhiS_m(2, x)
ctx_calc_m = ctx.calc_m(2, x)

# --- Compare results ---
def assert_close(name, a, b, tol=1e-14):
    if isinstance(a, (list, tuple)):
        for i, (ai, bi) in enumerate(zip(a, b)):
            assert_close(f"{name}[{i}]", ai, bi, tol)
        return
    if math.isnan(a) and math.isnan(b):
        return
    diff = abs(a - b)
    rel = diff / max(abs(a), 1e-300)
    status = "PASS" if rel < tol else "FAIL"
    if status == "FAIL":
        print(f"  {status}: {name}: global={a}, ctx={b}, rel_diff={rel}")
    assert rel < tol, f"{name}: global={a}, ctx={b}, rel_diff={rel}"

print("Comparing calc_PhiS...")
assert_close("PhiS", global_PhiS, ctx_PhiS)

print("Comparing calc...")
assert_close("calc.PhiS", global_calc[0], ctx_calc[0])
assert_close("calc.dPhiS", global_calc[1], ctx_calc[1])
assert_close("calc.d2PhiS", global_calc[2], ctx_calc[2])
assert_close("calc.src", global_calc[3], ctx_calc[3])

print("Comparing calc_PhiS_m...")
assert_close("PhiS_m", global_PhiS_m, ctx_PhiS_m)

print("Comparing calc_m...")
assert_close("calc_m.PhiS", global_calc_m[0], ctx_calc_m[0])
assert_close("calc_m.dPhiS", global_calc_m[1], ctx_calc_m[1])
assert_close("calc_m.d2PhiS", global_calc_m[2], ctx_calc_m[2])
assert_close("calc_m.src", global_calc_m[3], ctx_calc_m[3])

# --- Test independence of two contexts ---
print("\nTesting context independence...")
ctx1 = es.EffsourceEquatorialContext(M, a)
ctx2 = es.EffsourceEquatorialContext(M, 0.9)  # different spin

ctx1.set_particle(x_p, E, L, ur_p)

# Different particle params for ctx2
E2 = 0.97
L2 = 3.8
x_p2 = es.make_coordinate(t=0.0, r=8.0, theta=theta_p, phi=0.5)
ctx2.set_particle(x_p2, E2, L2, 0.05)

# Compute on ctx2 (should NOT affect ctx1)
x2 = es.make_coordinate(t=0.0, r=8.1, theta=theta_p + 0.01, phi=0.52)
_ = ctx2.calc(x2)

# ctx1 should still give the same result as the global API
ctx1_calc = ctx1.calc(x)
assert_close("independence.PhiS", global_calc[0], ctx1_calc[0])
assert_close("independence.dPhiS", global_calc[1], ctx1_calc[1])
assert_close("independence.d2PhiS", global_calc[2], ctx1_calc[2])
assert_close("independence.src", global_calc[3], ctx1_calc[3])

print("\nAll tests passed!")
