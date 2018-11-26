#ifndef __COHERENT_H__
#define __COHERENT_H__

#include <complex.h>

#include "bigpsi.h"
#include "config.h"

double complex Gminus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Gplus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Gzero(struct bigpsi_t *psi,int L,struct configuration_t *config);

double complex Lambdaplus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Lambdaminus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Lambdazero(struct bigpsi_t *psi,int L,struct configuration_t *config);

double complex Vbeta(struct bigpsi_t *psi,int L,struct configuration_t *config);

int sc_time_evolution_coherent(double t,const double y[],double dydt[],void *p);
int laser_time_evolution_coherent(double t,const double y[],double dydt[],void *data);

double norm_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config);
double norm_qp_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config);
double norm_phonons_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config);

double complex fcosthetasquared_coherent(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config);
double complex costhetasquared_coherent(struct bigpsi_t *psi,struct configuration_t *config);
double complex costheta2d_coherent(struct bigpsi_t *psi,struct configuration_t *config);
double complex overlapS_coherent(struct bigpsi_t *psi,double *y0,double t0,struct configuration_t *config);

#endif //__COHERENT_H__
