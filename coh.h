#ifndef __COH_H__
#define __COH_H__

#include "bigpsi.h"
#include "config.h"

double complex Gminus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Gplus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Gzero(struct bigpsi_t *psi,int L,struct configuration_t *config);

double complex Lambdaplus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Lambdaminus(struct bigpsi_t *psi,int L,struct configuration_t *config);
double complex Lambdazero(struct bigpsi_t *psi,int L,struct configuration_t *config);

double complex Vbeta(struct bigpsi_t *psi,int L,struct configuration_t *config);

int sc_time_evolution_coh(double t,const double y[],double dydt[],void *p);
int big_sc_time_evolution_coh(double t,const double y[],double dydt[],void *data);

#endif //__COH_H__
