/* Helpers shared between kerr-equatorial-ctx.c and kerr-equatorial-split-ctx.c. */
#ifndef EFFSOURCE_EQUATORIAL_INTERNAL_H
#define EFFSOURCE_EQUATORIAL_INTERNAL_H

#include <complex.h>

extern const double effsource_ReEI[21][2][5][27];
extern const double effsource_ImEI[21][2][5][27];

void effsource_ellPQ_derivs_from_C1l(long double C1, long double *out);

long double _Complex effsource_box_m_channel(long double a, long double r,
    long double theta,
    long double _Complex dPr, long double _Complex dPth,
    long double _Complex dtt, long double _Complex dtph, long double _Complex drr,
    long double _Complex dthth, long double _Complex dphph);

#endif
