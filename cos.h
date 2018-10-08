#ifndef __COS_H__
#define __COS_H__

#include "bigpsi.h"
#include "config.h"

double complex Across(struct bigpsi_t *psi,int L,int Lprime,int n,struct configuration_t *config);
double complex costheta2d(struct bigpsi_t *psi,struct configuration_t *config);
double complex costhetasquared(struct bigpsi_t *psi,struct configuration_t *config);
double complex costhetasquaredLLprime(struct bigpsi_t *psi,int L,int Lprime,struct configuration_t *config);

#endif //__COS_H__
