%module effsource_equatorial

%{
#include "effsource_equatorial.h"
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

/* Wrap the equatorial context struct (opaque fields — users access via API functions) */
struct effsource_equatorial_ctx {
    struct coordinate xp;
    double M, a;
    double A0060, A0061, A0080, A0081, A0240, A0241, A0260, A0261, A0420, A0421,
           A0440, A0441, A0600, A0601, A0620, A0621, A0800, A0801, A1060, A1061,
           A1080, A1240, A1241, A1260, A1420, A1421, A1440, A1600, A1601, A1620,
           A1800, A2040, A2041, A2060, A2061, A2220, A2221, A2240, A2241, A2400,
           A2401, A2420, A2421, A2600, A2601, A3040, A3041, A3060, A3220, A3221,
           A3240, A3400, A3401, A3420, A3600, A4020, A4021, A4040, A4041, A4200,
           A4201, A4220, A4221, A4400, A4401, A5020, A5021, A5040, A5200, A5201,
           A5220, A5400, A6000, A6001, A6020, A6021, A6200, A6201, A7000, A7001,
           A7020, A7200, A8000, A8001, A9000;
    double alpha20, alpha02, beta, c;
    double rt, urt, rtt, urtt, phit, phitt;
    double dAdt0060, dAdt0061, dAdt0080, dAdt0081, dAdt0240, dAdt0241, dAdt0260,
           dAdt0261, dAdt0420, dAdt0421, dAdt0440, dAdt0441, dAdt0600, dAdt0601,
           dAdt0620, dAdt0621, dAdt0800, dAdt0801, dAdt1060, dAdt1061, dAdt1080,
           dAdt1240, dAdt1241, dAdt1260, dAdt1420, dAdt1421, dAdt1440, dAdt1600,
           dAdt1601, dAdt1620, dAdt1800, dAdt2040, dAdt2041, dAdt2060, dAdt2061,
           dAdt2220, dAdt2221, dAdt2240, dAdt2241, dAdt2400, dAdt2401, dAdt2420,
           dAdt2421, dAdt2600, dAdt2601, dAdt3040, dAdt3041, dAdt3060, dAdt3220,
           dAdt3221, dAdt3240, dAdt3400, dAdt3401, dAdt3420, dAdt3600, dAdt4020,
           dAdt4021, dAdt4040, dAdt4041, dAdt4200, dAdt4201, dAdt4220, dAdt4221,
           dAdt4400, dAdt4401, dAdt5020, dAdt5021, dAdt5040, dAdt5200, dAdt5201,
           dAdt5220, dAdt5400, dAdt6000, dAdt6001, dAdt6020, dAdt6021, dAdt6200,
           dAdt6201, dAdt7000, dAdt7001, dAdt7020, dAdt7200, dAdt8000, dAdt8001,
           dAdt9000;
    double dalphadt20, dalphadt02, dbetadt, dcdt;
    double dC1_dt02, dC1_dt10, dC1_dt20;
    double d2Adt20060, d2Adt20061, d2Adt20080, d2Adt20081, d2Adt20240, d2Adt20241,
           d2Adt20260, d2Adt20261, d2Adt20420, d2Adt20421, d2Adt20440, d2Adt20441,
           d2Adt20600, d2Adt20601, d2Adt20620, d2Adt20621, d2Adt20800, d2Adt20801,
           d2Adt21060, d2Adt21061, d2Adt21080, d2Adt21240, d2Adt21241, d2Adt21260,
           d2Adt21420, d2Adt21421, d2Adt21440, d2Adt21600, d2Adt21601, d2Adt21620,
           d2Adt21800, d2Adt22040, d2Adt22041, d2Adt22060, d2Adt22061, d2Adt22220,
           d2Adt22221, d2Adt22240, d2Adt22241, d2Adt22400, d2Adt22401, d2Adt22420,
           d2Adt22421, d2Adt22600, d2Adt22601, d2Adt23040, d2Adt23041, d2Adt23060,
           d2Adt23220, d2Adt23221, d2Adt23240, d2Adt23400, d2Adt23401, d2Adt23420,
           d2Adt23600, d2Adt24020, d2Adt24021, d2Adt24040, d2Adt24041, d2Adt24200,
           d2Adt24201, d2Adt24220, d2Adt24221, d2Adt24400, d2Adt24401, d2Adt25020,
           d2Adt25021, d2Adt25040, d2Adt25200, d2Adt25201, d2Adt25220, d2Adt25400,
           d2Adt26000, d2Adt26001, d2Adt26020, d2Adt26021, d2Adt26200, d2Adt26201,
           d2Adt27000, d2Adt27001, d2Adt27020, d2Adt27200, d2Adt28000, d2Adt28001,
           d2Adt29000;
    double d2alphadt220, d2alphadt202, d2betadt2, d2cdt2;
    double d2C1_dt200, d2C1_dt202, d2C1_dt210, d2C1_dt220;
};

// /* Simple functions - wrap directly (global API) */
// void effsource_init(double M, double a);
// void effsource_set_particle(struct coordinate *x_p, double e, double l, double ur_p);

// /* Expose the raw output-pointer functions for use by the Python wrappers below */
// void effsource_PhiS(struct coordinate *x, double *PhiS);
// void effsource_calc(struct coordinate *x,
//     double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);
// void effsource_PhiS_m(int m, struct coordinate *x, double *PhiS);
// void effsource_calc_m(int m, struct coordinate *x,
//     double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);

/* Context-based API */
struct effsource_equatorial_ctx * effsource_equatorial_create(double M, double a);
void effsource_equatorial_free(struct effsource_equatorial_ctx * ctx);
void effsource_equatorial_ctx_set_particle(struct effsource_equatorial_ctx * ctx, struct coordinate *x_p, double E, double L, double ur);
void effsource_equatorial_ctx_PhiS(struct effsource_equatorial_ctx * ctx, struct coordinate *x, double *PhiS);
void effsource_equatorial_ctx_calc(struct effsource_equatorial_ctx * ctx, struct coordinate *x,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);
void effsource_equatorial_ctx_PhiS_m(struct effsource_equatorial_ctx * ctx, int m, struct coordinate *x, double *PhiS);
void effsource_equatorial_ctx_calc_m(struct effsource_equatorial_ctx * ctx, int m, struct coordinate *x,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);

/* Offset-based entry points (avoid near-particle coordinate cancellation) */
void effsource_equatorial_ctx_PhiS_offset(struct effsource_equatorial_ctx * ctx,
    double dr, double dtheta, double dphi, double *PhiS);
void effsource_equatorial_ctx_calc_offset(struct effsource_equatorial_ctx * ctx,
    double dr, double dtheta, double dphi,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);
void effsource_equatorial_ctx_PhiS_m_offset(struct effsource_equatorial_ctx * ctx,
    int m, double dr, double dtheta, double *PhiS);
void effsource_equatorial_ctx_calc_m_offset(struct effsource_equatorial_ctx * ctx,
    int m, double dr, double dtheta,
    double *PhiS, double *dPhiS_dx, double *d2PhiS_dx2, double *src);

/* Kernel split: PhiS_m = G0 + GL*ln(alpha), out = {ReG0, ImG0, ReGL, ImGL} */
void effsource_equatorial_ctx_PhiS_m_split(struct effsource_equatorial_ctx * ctx,
    int m, double dr, double dtheta, double *out);

/* Kernel quadratic-form coefficients: out = {alpha20, alpha02, beta, c} */
void effsource_equatorial_ctx_get_alpha(struct effsource_equatorial_ctx * ctx,
    double *out);

/* Kernel split of calc_m: 4 channels {A, L, P1, P2} per component,
   value = A + L*ln(alpha) + P1/alpha + P2/alpha^2 */
void effsource_equatorial_ctx_calc_m_split(struct effsource_equatorial_ctx * ctx,
    int m, double dr, double dtheta,
    double *PhiS_s, double *dPhiS_s, double *d2PhiS_s, double *src_s);

/* Utility to disable GSL error handler (prevents aborts on roundoff errors) */
%inline %{
void disable_gsl_error_handler(void) {
    gsl_set_error_handler_off();
}
%}

/* Python convenience wrappers */
%pythoncode %{
def make_coordinate(t=0.0, r=10.0, theta=0.0, phi=0.0):
    """Create a coordinate struct with the given values."""
    x = coordinate()
    x.t = t
    x.r = r
    x.theta = theta
    x.phi = phi
    return x

class EffsourceEquatorialContext:
    """Thread-safe context for equatorial effective source computations.

    Each instance holds its own state (particle position, coefficients, etc.),
    so multiple instances can be used in parallel without interference.
    """

    def __init__(self, M, a):
        self._ctx = effsource_equatorial_create(M, a)

    def __del__(self):
        if hasattr(self, '_ctx') and self._ctx is not None:
            effsource_equatorial_free(self._ctx)
            self._ctx = None

    def set_particle(self, x_p, E, L, ur):
        effsource_equatorial_ctx_set_particle(self._ctx, x_p, E, L, ur)

    def calc_PhiS(self, x):
        """Compute singular field at point x.

        Returns: float
        """
        buf = doubleArray(1)
        effsource_equatorial_ctx_PhiS(self._ctx, x, buf.cast())
        return buf[0]

    def calc_PhiS_m(self, m, x):
        """Compute m-mode of singular field at point x.

        Returns: (Re, Im) tuple
        """
        buf = doubleArray(2)
        effsource_equatorial_ctx_PhiS_m(self._ctx, m, x, buf.cast())
        return buf[0], buf[1]

    def calc(self, x):
        """Compute effective source at point x.

        Returns: (PhiS, dPhiS_dx, d2PhiS_dx2, src)
        """
        _PhiS = doubleArray(1)
        _dPhiS = doubleArray(4)
        _d2PhiS = doubleArray(10)
        _src = doubleArray(1)
        effsource_equatorial_ctx_calc(self._ctx, x, _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
        return (
            _PhiS[0],
            [_dPhiS[i] for i in range(4)],
            [_d2PhiS[i] for i in range(10)],
            _src[0],
        )

    def calc_m(self, m, x):
        """Compute m-mode effective source at point x.

        Returns: (PhiS, dPhiS, d2PhiS, src)
        """
        _PhiS = doubleArray(2)
        _dPhiS = doubleArray(8)
        _d2PhiS = doubleArray(20)
        _src = doubleArray(2)
        effsource_equatorial_ctx_calc_m(self._ctx, m, x, _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
        return (
            [_PhiS[i] for i in range(2)],
            [_dPhiS[i] for i in range(8)],
            [_d2PhiS[i] for i in range(20)],
            [_src[i] for i in range(2)],
        )

    def calc_PhiS_offset(self, dr, dtheta, dphi=0.0):
        """Singular field at an offset (dr, dtheta, dphi) from the particle."""
        buf = doubleArray(1)
        effsource_equatorial_ctx_PhiS_offset(self._ctx, dr, dtheta, dphi, buf.cast())
        return buf[0]

    def calc_PhiS_m_offset(self, m, dr, dtheta):
        """m-mode singular field at an offset (dr, dtheta) from the particle.

        Returns: (Re, Im) tuple
        """
        buf = doubleArray(2)
        effsource_equatorial_ctx_PhiS_m_offset(self._ctx, m, dr, dtheta, buf.cast())
        return buf[0], buf[1]

    def calc_offset(self, dr, dtheta, dphi=0.0):
        """Effective source at an offset (dr, dtheta, dphi) from the particle.

        Returns: (PhiS, dPhiS_dx, d2PhiS_dx2, src)
        """
        _PhiS = doubleArray(1)
        _dPhiS = doubleArray(4)
        _d2PhiS = doubleArray(10)
        _src = doubleArray(1)
        effsource_equatorial_ctx_calc_offset(self._ctx, dr, dtheta, dphi,
                                             _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
        return (
            _PhiS[0],
            [_dPhiS[i] for i in range(4)],
            [_d2PhiS[i] for i in range(10)],
            _src[0],
        )

    def calc_PhiS_m_split(self, m, dr, dtheta):
        """Kernel split of the m-mode singular field at an offset:
        PhiS_m = G0 + GL*ln(alpha), alpha = alpha20*dr^2 + alpha02*dtheta^2,
        with G0, GL analytic across the particle.

        Returns: ((ReG0, ImG0), (ReGL, ImGL))
        """
        buf = doubleArray(4)
        effsource_equatorial_ctx_PhiS_m_split(self._ctx, m, dr, dtheta, buf.cast())
        return (buf[0], buf[1]), (buf[2], buf[3])

    def get_alpha(self):
        """Kernel quadratic-form coefficients for the current particle.

        Returns: (alpha20, alpha02, beta, c)
        """
        buf = doubleArray(4)
        effsource_equatorial_ctx_get_alpha(self._ctx, buf.cast())
        return buf[0], buf[1], buf[2], buf[3]

    def calc_m_split(self, m, dr, dtheta):
        """Kernel split of calc_m at an offset (dr, dtheta): every component
        of (PhiS, dPhiS, d2PhiS, src) as seven channels {A, L, P1, P2, P3, P4,
        P5} with value = A + L*ln(alpha) + P1/alpha + ... + P5/alpha^5 and all
        channels analytic across the particle (the P-side hidden poles are
        extracted in closed form; field 1/alpha^3, 1st derivs 1/alpha^4, 2nd
        derivs / src 1/alpha^5).

        Returns: (PhiS_s, dPhiS_s, d2PhiS_s, src_s) where each entry is a
        list of 7 channel blocks with calc_m_offset's component layout
        (mixed second derivatives that calc_m leaves NAN are 0 here).
        """
        _PhiS = doubleArray(14)
        _dPhiS = doubleArray(56)
        _d2PhiS = doubleArray(140)
        _src = doubleArray(14)
        effsource_equatorial_ctx_calc_m_split(self._ctx, m, dr, dtheta,
                                              _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
        return (
            [[_PhiS[2 * ch + i] for i in range(2)] for ch in range(7)],
            [[_dPhiS[8 * ch + i] for i in range(8)] for ch in range(7)],
            [[_d2PhiS[20 * ch + i] for i in range(20)] for ch in range(7)],
            [[_src[2 * ch + i] for i in range(2)] for ch in range(7)],
        )

    def calc_m_offset(self, m, dr, dtheta):
        """m-mode effective source at an offset (dr, dtheta) from the particle.

        Returns: (PhiS, dPhiS, d2PhiS, src)
        """
        _PhiS = doubleArray(2)
        _dPhiS = doubleArray(8)
        _d2PhiS = doubleArray(20)
        _src = doubleArray(2)
        effsource_equatorial_ctx_calc_m_offset(self._ctx, m, dr, dtheta,
                                               _PhiS.cast(), _dPhiS.cast(), _d2PhiS.cast(), _src.cast())
        return (
            [_PhiS[i] for i in range(2)],
            [_dPhiS[i] for i in range(8)],
            [_d2PhiS[i] for i in range(20)],
            [_src[i] for i in range(2)],
        )
%}
