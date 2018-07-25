#ifndef __OBSERVABLES_H__
#define __OBSERVABLES_H__

#include <stdio.h>
#include <math.h>

#include "bigpsi.h"
#include "config.h"

struct dint_container_t
{
	struct interpolation_t *intre1;
	struct interpolation_t *intim1;

	struct interpolation_t *intre2;
	struct interpolation_t *intim2;

#define DINT_MODE_PLAIN		(21)
#define DINT_MODE_OMEGAK	(22)
#define DINT_MODE_VK		(23)
#define DINT_MODE_VK0		(24)

	int L,Lprime,mode;
	
	struct configuration_t *config;
	
	double localdensity;
};

double complex Dcross(struct bigpsi_t *psi,int L,int Lprime,int n,int nprime,int mode,struct configuration_t *config);
double complex Dsingle(struct bigpsi_t *psi,int L,int Lprime,int n,int mode,struct configuration_t *config);
double complex Ecross(struct bigpsi_t *psi,int L,int Lprime,int n,double *y0,double t0,struct configuration_t *config);

double eta_sigma(int L, int lambda,int n,int nu);

double complex molecular_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config);
double complex bosons_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config);
double complex total_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config);

double complex rcr(int L,struct bigpsi_t *psi,struct configuration_t *config);
double complex overlapS(struct bigpsi_t *psi,double *y0,double t0,struct configuration_t *config);
double torque(struct bigpsi_t *psi,int L,int M,struct configuration_t *config);

#endif //__OBSERVABLES_H__
