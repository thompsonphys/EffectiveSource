%module effsource_circular

%{
#include "effsource.h"
#include <gsl/gsl_errno.h>
%}

%include "carrays.i"
%array_class(double, doubleArray);

/* Wrap the coordinate struct */
struct coordinate {
    double r;
    double theta;
    double phi;
    double t;
};

/* Simple functions - wrap directly */
void effsource_init(double M, double a);
void effsource_set_particle(struct coordinate *x_p, double e, double l, double ur_p);

/* Expose the raw output-pointer functions for use by the Python wrappers below */
void effsource_PhiS(struct coordinate *x, double *PhiS);
void effsource_calc(struct coordinate *x,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);
void effsource_PhiS_m(int m, struct coordinate *x, double *PhiS);
void effsource_calc_m(int m, struct coordinate *x,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);

/* Utility to disable GSL error handler (prevents aborts on roundoff errors) */
%inline %{
void disable_gsl_error_handler(void) {
    gsl_set_error_handler_off();
}
%}

/* Python convenience wrappers */
%pythoncode %{
def make_coordinate(t=0.0, r=0.0, theta=0.0, phi=0.0):
    """Create a coordinate struct with the given values."""
    x = coordinate()
    x.t = t
    x.r = r
    x.theta = theta
    x.phi = phi
    return x

def calc_PhiS(x):
    """Compute singular field at point x.

    Returns: float
    """
    buf = doubleArray(1)
    effsource_PhiS(x, buf.cast())
    return buf[0]

def calc_PhiS_m(m, x):
    """Compute m-mode of singular field at point x.

    Returns: (Re, Im) tuple
    """
    buf = doubleArray(2)
    effsource_PhiS_m(m, x, buf.cast())
    return buf[0], buf[1]

def calc(x):
    """Compute effective source at point x.

    Returns: (PhiS, dPhiS_dx, d2PhiS_dx2, src)
        PhiS: float - singular field
        dPhiS_dx: list[4] - first derivatives [t, r, theta, phi]
        d2PhiS_dx2: list[10] - second derivatives
        src: float - effective source
    """
    _PhiS = doubleArray(1)
    _dPhiS = doubleArray(4)
    _d2PhiS = doubleArray(10)
    _src = doubleArray(1)
    effsource_calc(x, _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
    return (
        _PhiS[0],
        [_dPhiS[i] for i in range(4)],
        [_d2PhiS[i] for i in range(10)],
        _src[0],
    )

def calc_m(m, x):
    """Compute m-mode effective source at point x.

    Returns: (PhiS, dPhiS, d2PhiS, src)
        PhiS: list[2] - [Re, Im] singular field
        dPhiS: list[8] - first derivatives (Re/Im interleaved)
        d2PhiS: list[20] - second derivatives (Re/Im interleaved)
        src: list[2] - [Re, Im] effective source
    """
    _PhiS = doubleArray(2)
    _dPhiS = doubleArray(8)
    _d2PhiS = doubleArray(20)
    _src = doubleArray(2)
    effsource_calc_m(m, x, _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
    return (
        [_PhiS[i] for i in range(2)],
        [_dPhiS[i] for i in range(8)],
        [_d2PhiS[i] for i in range(20)],
        [_src[i] for i in range(2)],
    )
%}
