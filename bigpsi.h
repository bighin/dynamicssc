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
	
	/*
		Saved norms for normalization
	*/
	
	bool have_normalization_snapshot;
	double *normalization_snapshot;
};

double Q(int L,int Lprime,int M,int N);
bool laser_is_on(struct configuration_t *config,double t);

int big_sc_time_evolution(double t,const double y[],double dydt[],void *data);

#define BIGPSI_INIT_FROM_CONFIG	(31)
#define BIGPSI_INIT_FROM_VALUES	(32)

struct bigpsi_t *bigpsi_init(struct configuration_t *config,int mode,int L,int M);
void bigpsi_fini(struct bigpsi_t *psi);

void bigpsi_serialize(struct bigpsi_t *psi,FILE *out);
struct bigpsi_t *bigpsi_deserialize(FILE *in,struct configuration_t *config);

double get_aos(struct bigpsi_t *psi);
double total_norm(struct bigpsi_t *psi);
double total_norm_qp(struct bigpsi_t *psi);
double total_norm_phonons(struct bigpsi_t *psi);

void bigpsi_normalize(struct bigpsi_t *psi,double *previousnorm);
void bigpsi_apply_step(struct bigpsi_t *psi,double ti,double *previousnorm,struct configuration_t *config);

#endif //__BIGPSI_H__
