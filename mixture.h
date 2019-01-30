#ifndef __MIXTURE_H__
#define __MIXTURE_H__

#include <complex.h>

#include "config.h"

struct info_t
{
	double t;
	double intensity,bath_intensity;
	double totalnorm,totalnorm_qp,totalnorm_ph;

	double *norms;
	double *norms_qp;
	double *norms_ph;

	double complex cos2d,cossquared;
	double aos;
	double torque;

	double j2,l2,lambda2;
};

int do_run(int L,int M,struct info_t *info,bool silent,struct configuration_t *config);
void do_mixture(struct configuration_t *config);

#endif //__MIXTURE_H__
