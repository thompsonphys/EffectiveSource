/*******************************************************************************
 * Copyright (C) 2011 Barry Wardell
 ******************************************************************************/

struct coordinate {
  double r;
  double theta;
  double phi;
  double t;
};

double xFunc(double a, double p, double e);

void effsource_init(double M, double a);
void effsource_set_particle(struct coordinate * x_p, double E, double L, double ur_p);

void effsource_PhiS(struct coordinate * x, double * PhiS);
void effsource_calc(struct coordinate * x,
  double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);

void effsource_PhiS_m(int m, struct coordinate * x, double * PhiS);
void effsource_calc_m(int m, struct coordinate * x,
  double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);

/* Thread-safe context-based API */
struct effsource_ctx {
  struct coordinate xp;
  double M, a;

  /* Series expansion coefficients (set by effsource_ctx_set_particle) */
  double A006, A008, A024, A026, A042, A044, A060, A062, A080,
         A106, A108, A124, A126, A142, A144, A160, A162, A180,
         A204, A206, A222, A224, A240, A242, A260,
         A304, A306, A322, A324, A340, A342, A360,
         A402, A404, A420, A422, A440,
         A502, A504, A520, A522, A540,
         A600, A602, A620,
         A700, A702, A720,
         A800, A900;
  double alpha20, alpha02, beta;
};

struct effsource_ctx * effsource_create(double M, double a);
void effsource_free(struct effsource_ctx * ctx);
void effsource_ctx_set_particle(struct effsource_ctx * ctx, struct coordinate * x_p, double E, double L, double ur_p);
void effsource_ctx_PhiS(struct effsource_ctx * ctx, struct coordinate * x, double * PhiS);
void effsource_ctx_calc(struct effsource_ctx * ctx, struct coordinate * x,
  double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);
void effsource_ctx_PhiS_m(struct effsource_ctx * ctx, int m, struct coordinate * x, double * PhiS);
void effsource_ctx_calc_m(struct effsource_ctx * ctx, int m, struct coordinate * x,
  double * PhiS, double * dPhiS_dx, double * d2PhiS_dx2, double * src);
