/* Shared declarations for the original global-scope API. Included by
   effsource_original.i and effsource_original_agm.i, which supply the
   %module line and differ only in which kerr-equatorial.c they link. */

%{
#include "effsource.h"
#include <gsl/gsl_errno.h>
%}

%include "carrays.i"
%array_class(double, doubleArray);

struct coordinate {
    double r;
    double theta;
    double phi;
    double t;
};

void effsource_init(double M, double a);
void effsource_set_particle(struct coordinate * x_p, double E, double L, double ur_p);
void effsource_PhiS(struct coordinate * x, double * PhiS);
void effsource_calc(struct coordinate * x,
    double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);
void effsource_PhiS_m(int m, struct coordinate * x, double * PhiS);
void effsource_calc_m(int m, struct coordinate * x,
    double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);

%inline %{
void disable_gsl_error_handler(void) {
    gsl_set_error_handler_off();
}
%}

%pythoncode %{
def make_coordinate(t=0.0, r=10.0, theta=0.0, phi=0.0):
    """(t, r, theta, phi) -> coordinate struct."""
    x = coordinate()
    x.t = t
    x.r = r
    x.theta = theta
    x.phi = phi
    return x
%}
