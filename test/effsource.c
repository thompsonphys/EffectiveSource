#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include "effsource.h"
#include "decompose.h"

/* Return the value of the singular field at the point x */
double PhiS_calc(struct coordinate * x)
{
  double PhiS;
  effsource_PhiS(x, &PhiS);
  return PhiS;
}

/* Return the value of the effective source at the point x */
double src_calc(struct coordinate * x)
{
  double PhiS, dPhiS_dx[4], d2PhiS_dx2[10], src;
  effsource_calc(x, &PhiS, dPhiS_dx, d2PhiS_dx2, &src);
  return src;
}

double xFunc(double a, double p, double e)
{
  double x, F, N, C, signA;
  double p2, p3, e2p3;

  signA = a > 0.0 ? 1.0 : -1.0;

  p2 = p * p;
  p3 = p2 * p;
  e2p3 = 3.0 + e*e;

  C = (a * a - p);
  C *= C;
  N = 2.0 / p * (-p2 + (e2p3 - a * a) * p - a * a * (1 + 3.0 * e * e) );
  F = 1 / p3 * (
    p3 - 2.0*e2p3*p2 + e2p3*e2p3*p - 4.0 * a * a * (1.0 - e * e) * (1.0 - e * e)
  );

  x = sqrt((-N-signA*sqrt(N*N-4.0*F*C))/(2.0*F));

  return x;
}

int main(int argc, char* argv[])
{
  FILE *inputFile, *outputFile;
  char buf[0x100];

  /* set m-mode to compute */
  int m = 5;

  /* Mass and spin of the central black hole */
  const double a = 0.5;
  const double M = 1.0;
  const double p = 9.9;
  const double e = 0.1;
  effsource_init(M, a);

  /* Orbital parameters */
  double r_p  = 10.0;
  double t = 0.0;
  double phi_p = 0.0;
  double E, L, X, ur;
  struct coordinate xp;
  struct coordinate x;

  /* set energy and angular momentum */
  X = xFunc(a,p,e);
  E = sqrt(
    1 - (M / p) * (1.0 - e * e) * (1 - X * X / (p * p) * (1.0 - e * e) )
  );
  L = X + a * E;

  /* Output the singular field for the m=2 mode in the r-theta plane */
  double r = 10.0;
  double theta = M_PI_2 - 0.01;
  double PhiS[2], dPhiS[8], ddPhiS[20], src[2];

  x.r = r;
  x.theta = theta;

  // inputFile = fopen("/Users/jthompson/projects/sf/effsource_fd/trajectory_data_t.dat","r");
  // if(inputFile == NULL){
  //   printf("No file found.");
  // }

  // snprintf(buf, sizeof(buf), "output_data/output_data_m%d.dat", m);
  // outputFile = fopen(buf,"w");
  // if(outputFile == NULL){
  //   printf("No file to write to.");
  // }

  // fprintf(outputFile, "# m = %d\n",m);
  // fprintf(outputFile, "# a = %lf, p = %lf, e = %lf, x = %lf\n",a,p,e,1.0);
  // fprintf(outputFile, "# t\tr\ttheta\tphi\tRePhi\tImPhi\n");

  // if(inputFile != NULL){
  //   while(fscanf(inputFile,"%lf\t%lf\t%lf\t%lf", &t, &r_p, &phi_p,&ur) != EOF){

  //     xp.t = t;
  //     xp.r = r_p;
  //     xp.theta = M_PI_2;
  //     xp.phi = phi_p;

  //     effsource_set_particle(&xp, E, L, ur);

  //     /* The point where we measure the singular field/effective source */
  //     struct coordinate x = {0.0, r, theta, 0.0};

  //     /* Disable the GSL error handler so that it doesn't abort due to roundoff errors */
  //     gsl_set_error_handler_off();

  //     effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);

  //     fprintf(outputFile, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
  //       t, x.r, x.theta, x.phi, src[0], src[1]);
  //   }
  // }

  inputFile = fopen("/Users/jthompson/projects/sf/effsource_fd/trajectory_data_chi.dat","r");
  if(inputFile == NULL){
    printf("No file found.");
  }

  snprintf(buf, sizeof(buf), "output_data/output_data_chi_m%d.dat", m);
  outputFile = fopen(buf,"w");
  if(outputFile == NULL){
    printf("No file to write to.");
  }

  fprintf(outputFile, "# m = %d\n",m);
  fprintf(outputFile, "# a = %lf, p = %lf, e = %lf, x = %lf\n",a,p,e,1.0);
  fprintf(outputFile, "# chi\tr\ttheta\tphi\tRePhi\tImPhi\n");

  if(inputFile != NULL){
    while(fscanf(inputFile,"%lf\t%lf\t%lf\t%lf", &t, &r_p, &phi_p,&ur) != EOF){

      xp.t = 0.0;
      xp.r = r_p;
      xp.theta = M_PI_2;
      xp.phi = phi_p;

      effsource_set_particle(&xp, E, L, ur);

      /* The point where we measure the singular field/effective source */
      struct coordinate x = {0.0, r, theta, 0.0};

      /* Disable the GSL error handler so that it doesn't abort due to roundoff errors */
      gsl_set_error_handler_off();

      effsource_calc_m(m, &x, PhiS, dPhiS, ddPhiS, src);

      fprintf(outputFile, "%.15g\t%.15g\t%.15g\t%.15g\t%.15g\t%.15g\n",
        t, x.r, x.theta, x.phi, src[0], src[1]);
    }
  }

  fclose(inputFile);
  fclose(outputFile);

  return 0;
}