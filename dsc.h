#ifndef __DSC_H__
#define __DSC_H__

#include <complex.h>
#include <gsl/gsl_odeiv2.h>

#include "config.h"

/*
	Linux doesn't have this definition for some reason
*/

#ifndef M_PI
#define M_PI	(3.14159265358979323846f)
#endif


double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3);
double cg(int j1,int m1,int j2,int m2,int j,int m);

double omegak(double k);
double U0(double n,double k,struct configuration_t *config);
double U2(double n,double k,struct configuration_t *config);
double V0(double k,double n,struct configuration_t *config);
double V2(double k,double n,struct configuration_t *config);

double complex timephase(double phase,double t,struct configuration_t *config);
double adiabatic_ramp(double t,struct configuration_t *config);
double W(double k,struct configuration_t *config);

struct params_t
{
	int L,M;

	gsl_odeiv2_driver *driver;
	gsl_odeiv2_system sys;

	double t;
	
	struct configuration_t *config;
};

struct container_t
{
	struct interpolation_t *intrexi2m2;
	struct interpolation_t *intimxi2m2;

	struct interpolation_t *intrexi2m1;
	struct interpolation_t *intimxi2m1;

	struct interpolation_t *intrexi20;
	struct interpolation_t *intimxi20;

	struct interpolation_t *intrexi21;
	struct interpolation_t *intimxi21;

	struct interpolation_t *intrexi22;
	struct interpolation_t *intimxi22;

	double t,localdensity;

	struct params_t *params;
};

double norm_qp(double t,const double y[],struct configuration_t *config);
double norm_phonons(double t,const double y[],struct configuration_t *config);
double norm(double t,const double y[],struct configuration_t *config);

double complex Aplus(double t,const double y[],struct params_t *params,double localdensity);
double complex Aminus(double t,const double y[],struct params_t *params,double localdensity);

int sc_time_evolution(double t,const double y[],double dydt[],void *p);

#endif //__DSC_H__
