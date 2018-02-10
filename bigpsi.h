#ifndef __BIGPSI_H__
#define __BIGPSI_H__

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "dsc.h"
#include "config.h"

struct bigpsi_t
{
	/*
		A pointer to the global configuration, as read from an INI file
	*/

	struct configuration_t *config;

	/*
		A dynamically-allocated array of param_t structures, one for each L-state.
	*/

	struct params_t *params;

	/*
		The wavefunction and the time
	*/

	double *y;
	double t;

	int nrpsis,dim;

	/*
		The GSL integrator
	*/

	gsl_odeiv2_driver *driver;
	gsl_odeiv2_system sys;
};

int big_sc_time_evolution(double t,const double y[],double dydt[],void *data);

struct bigpsi_t *bigpsi_init(struct configuration_t *config,int L,int M);
void bigpsi_fini(struct bigpsi_t *psi);

double get_aos(struct bigpsi_t *psi);
double total_norm(struct bigpsi_t *psi);
double total_norm_qp(struct bigpsi_t *psi);
double total_norm_phonons(struct bigpsi_t *psi);

void bigpsi_normalize(struct bigpsi_t *psi,double *normalization_error);
void bigpsi_apply_step(struct bigpsi_t *psi,double ti,double *normalization_error);

#endif //__BIGPSI_H__
