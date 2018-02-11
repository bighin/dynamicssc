#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <sys/time.h>

#include <gsl/gsl_sf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

#include "cubature/cubature.h"
#include "libprogressbar/progressbar.h"
#include "auxx.h"
#include "config.h"
#include "dsc.h"

/*
        Tuning parameters for numerical integration.
*/

size_t maxEval=1024*1024;
double relError=1e-5;

/*
	Wigner's 3j symbol from the GNU Scientific Library
*/

double wigner3j(int l1,int l2,int l3,int m1,int m2,int m3)
{
	return gsl_sf_coupling_3j(2*l1,2*l2,2*l3,2*m1,2*m2,2*m3);
}

double cg(int j1,int m1,int j2,int m2,int j,int m)
{
	double phase;

	if(((j1-j2+m)%2)==0)
		phase=1.0f;
	else
		phase=-1.0f;

	return phase*sqrt(2.0f*j+1.0f)*wigner3j(j1,j2,j,m1,m2,-m);
}

/*
	Dispersion relations for free particles.
*/

double ek(double k)
{
	return k*k/2.0f;
}

/*
	The spectrum is obtained from experimental data (see Igor's notebook),
	and fitted using a Pade' approximant of order [5/5]
	with constraints \omega(0)=0 and \omega(\infty)=350.
*/

double omegak(double k)
{
	double num,den;

	/*
		Padè approximant for the dispersion relation
	*/

	num=k*(26565.2+k*(-1497.82+k*(129.012+k*(-7.37462+k*0.128625))));
	den=1154.95+k*(-73.3346+k*(4.19356+k*(0.168494+k*(-0.0193761+k*0.000367499))));

	return num/den;
}

/*
	Potentials
*/

double U0(double n,double k,struct configuration_t *config)
{
	double a,b,c;

	a=sqrt((8.0f*n*k*k*ek(k))/(omegak(k)));
	b=exp(-0.5f*k*k*pow(config->r0,2.0));
	c=4*M_PI*pow(config->r0,-3.0f);

	return config->u0*a*b/c;
}

double U2(double n,double k,struct configuration_t *config)
{
	double a,b,c,d,z;

	z=k*config->r2;
	a=sqrt((8.0f*n*k*k*ek(k))/(5.0f*omegak(k)));
	b=-2.0f*exp(-0.5*z*z)*z*(3.0f+z*z);
	c=3.0f*sqrt(2.0f*M_PI)*erf(z/sqrt(2.0f));
	d=8.0f*pow(k,3.0f)*M_PI;
	
	return config->u2*a*(b+c)/d;
}

double V0(double k,double n,struct configuration_t *config)
{
	if(k==0.0f)
		return 0.0f;

	return U0(n,k,config)*sqrtf(1.0f/(4.0f*M_PI));
}

double V2(double k,double n,struct configuration_t *config)
{
	if(k==0.0f)
		return 0.0f;

	return U2(n,k,config)*sqrtf(5.0f/(4.0f*M_PI));
}

double complex timephase(double phase,double t,struct configuration_t *config)
{
	double mintime=config->starttime;
	
	return cexp(I*phase*(t-mintime));
}

/*
	A simple ramp, using the function (1/2) * (1 + tanh((t - tstart) / delta))

	The time for going from a value x to a value y is:

	delta * (arctanh(1-2*y) - arctanh(1-2*x))

	so that the time to go from 0.1 to 0.9 is:

	delta * (arctanh(0.8) - arctanh(-0.8)) = 2 * delta * arctanh(0.8) \approx delta * 2.19

	If I choose delta = 4, then the ramp up time from 0.1 to 0.9 is about 9 (in units of B),
	which should be enough to be considered adiabatic.
*/

double adiabatic_ramp(double t,struct configuration_t *config)
{
	double rampcenter=config->rampcenter;
	double delta=config->rampdelta;

	return 0.5f*(1.0f+tanh((t-rampcenter)/delta));
}

int fnorm(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the 'external' variables */

        struct container_t *container=(struct container_t *)(fdata);

        /* The integration variables */

        double k;
	double complex xi2m2,xi2m1,xi20,xi21,xi22;

	k=x[0];

	xi2m2=get_point(container->intrexi2m2,k)+I*get_point(container->intimxi2m2,k);
	xi2m1=get_point(container->intrexi2m1,k)+I*get_point(container->intimxi2m1,k);
	xi20=get_point(container->intrexi20,k)+I*get_point(container->intimxi20,k);
	xi21=get_point(container->intrexi21,k)+I*get_point(container->intimxi21,k);
	xi22=get_point(container->intrexi22,k)+I*get_point(container->intimxi22,k);

	fval[0]=conj(xi2m2)*xi2m2+conj(xi2m1)*xi2m1+conj(xi20)*xi20+conj(xi21)*xi21+conj(xi22)*xi22;

	return 0;
}

double norm_qp(double t,const double y[],struct configuration_t *config)
{
	double complex g=y[0]+I*y[1];
	
	return sqrt(conj(g)*g);
}

double norm_phonons(double t,const double y[],struct configuration_t *config)
{
	double *x,*yre2m2,*yim2m2,*yre2m1,*yim2m1,*yre20,*yim20,*yre21,*yim21,*yre22,*yim22;
	double xmin,xmax,err;

	struct container_t container;

	if(config->freeevolution)
		return 0.0f;

	x=malloc(sizeof(double)*config->gridpoints);
	yre2m2=malloc(sizeof(double)*config->gridpoints);
	yim2m2=malloc(sizeof(double)*config->gridpoints);
	yre2m1=malloc(sizeof(double)*config->gridpoints);
	yim2m1=malloc(sizeof(double)*config->gridpoints);
	yre20=malloc(sizeof(double)*config->gridpoints);
	yim20=malloc(sizeof(double)*config->gridpoints);
	yre21=malloc(sizeof(double)*config->gridpoints);
	yim21=malloc(sizeof(double)*config->gridpoints);
	yre22=malloc(sizeof(double)*config->gridpoints);
	yim22=malloc(sizeof(double)*config->gridpoints);

	double res=0.0f;

        for(int c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;

                x[c]=k;
                yre2m2[c]=y[2+10*c];
                yim2m2[c]=y[2+10*c+1];
                yre2m1[c]=y[2+10*c+2];
                yim2m1[c]=y[2+10*c+3];
                yre20[c]=y[2+10*c+4];
                yim20[c]=y[2+10*c+5];
                yre21[c]=y[2+10*c+6];
                yim21[c]=y[2+10*c+7];
                yre22[c]=y[2+10*c+8];
                yim22[c]=y[2+10*c+9];
	}

	container.intrexi2m2=init_interpolation(x,yre2m2,config->gridpoints);
	container.intimxi2m2=init_interpolation(x,yim2m2,config->gridpoints);
	container.intrexi2m1=init_interpolation(x,yre2m1,config->gridpoints);
	container.intimxi2m1=init_interpolation(x,yim2m1,config->gridpoints);
	container.intrexi20=init_interpolation(x,yre20,config->gridpoints);
	container.intimxi20=init_interpolation(x,yim20,config->gridpoints);
	container.intrexi21=init_interpolation(x,yre21,config->gridpoints);
	container.intimxi21=init_interpolation(x,yim21,config->gridpoints);
	container.intrexi22=init_interpolation(x,yre22,config->gridpoints);
	container.intimxi22=init_interpolation(x,yim22,config->gridpoints);

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(1,fnorm,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,&res,&err);

	fini_interpolation(container.intrexi2m2);
	fini_interpolation(container.intimxi2m2);
	fini_interpolation(container.intrexi2m1);
	fini_interpolation(container.intimxi2m1);
	fini_interpolation(container.intrexi20);
	fini_interpolation(container.intimxi20);
	fini_interpolation(container.intrexi21);
	fini_interpolation(container.intimxi21);
	fini_interpolation(container.intrexi22);
	fini_interpolation(container.intimxi22);

	if(x) free(x);
	if(yre2m2) free(yre2m2);
	if(yim2m2) free(yim2m2);
	if(yre2m1) free(yre2m1);
	if(yim2m1) free(yim2m1);
	if(yre20) free(yre20);
	if(yim20) free(yim20);
	if(yre21) free(yre21);
	if(yim21) free(yim21);
	if(yre22) free(yre22);
	if(yim22) free(yim22);

	return sqrt(res);
}

double norm(double t,const double y[],struct configuration_t *config)
{
	return sqrt(pow(norm_qp(t,y,config),2.0f)+pow(norm_phonons(t,y,config),2.0f));
}

int fAplus(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase21;
	double t,localdensity;
	int L;
	
	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;
	
	k=x[0];
	phase21=timephase((L*(L+1.0f)+omegak(k)-2.0f),t,config);
	f=V2(k,localdensity,config)/omegak(k)*phase21*(get_point(container->intrexi21,k)+I*get_point(container->intimxi21,k));

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex Aplus(double t,const double y[],struct params_t *params,double localdensity)
{
	struct container_t container;
	struct configuration_t *config=params->config;

	double *x,*yre,*yim;
	double xmin,xmax,res[2],err[2];
	int c;

	x=malloc(sizeof(double)*config->gridpoints);
	yre=malloc(sizeof(double)*config->gridpoints);
	yim=malloc(sizeof(double)*config->gridpoints);

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;

                x[c]=k;
                yre[c]=y[2+10*c+6];
                yim[c]=y[2+10*c+7];
	}

        container.intrexi21=init_interpolation(x,yre,config->gridpoints);
        container.intimxi21=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fAplus,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi21);
	fini_interpolation(container.intimxi21);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int fAminus(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase2m1;
	double t,localdensity;
	int L;
	
	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;
	
	k=x[0];
	phase2m1=timephase((L*(L+1.0f)+omegak(k)-2.0f),t,config);
	f=V2(k,localdensity,config)/omegak(k)*phase2m1*(get_point(container->intrexi2m1,k)+I*get_point(container->intimxi2m1,k));

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex Aminus(double t,const double y[],struct params_t *params,double localdensity)
{
	struct container_t container;
	struct configuration_t *config=params->config;

	double *x,*yre,*yim;
	double xmin,xmax,res[2],err[2];
	int c;

	x=malloc(sizeof(double)*config->gridpoints);
	yre=malloc(sizeof(double)*config->gridpoints);
	yim=malloc(sizeof(double)*config->gridpoints);

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;

                x[c]=k;
                yre[c]=y[2+10*c+2];
                yim[c]=y[2+10*c+3];
	}

        container.intrexi2m1=init_interpolation(x,yre,config->gridpoints);
        container.intimxi2m1=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fAminus,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi2m1);
	fini_interpolation(container.intimxi2m1);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int sc_time_evolution(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	double complex localAplus,localAminus;
	double complex g,dgdt;
	double localdensity;
	int c,L;

	L=params->L;

	g=y[0]+I*y[1];

	if(config->freeevolution==true)
	{
		dgdt=I*L*(L+1)*g;

		dydt[0]=creal(dgdt);
		dydt[1]=cimag(dgdt);

		return GSL_SUCCESS;
	}
	
	if(config->ramp==true)
		localdensity=config->density*adiabatic_ramp(t,config);
	else
		localdensity=config->density;

	localAplus=Aplus(t,y,p,localdensity);
	localAminus=Aminus(t,y,p,localdensity);

	dgdt=I*(L*(L+1))*g+I*sqrt(6*L*(L+1))*(localAplus+localAminus);

	dydt[0]=creal(dgdt);
	dydt[1]=cimag(dgdt);

	if(L<=0)
		return GSL_SUCCESS;

	for(c=0;c<config->gridpoints;c++)
	{
		double gridstep=config->cutoff/config->gridpoints;
		double k=gridstep*c;

		double complex xi2m2,xi2m1,xi20,xi21,xi22;
		double complex dxi2m2dt,dxi2m1dt,dxi20dt,dxi21dt,dxi22dt;

		xi2m2=y[2+10*c]+I*y[2+10*c+1];
		xi2m1=y[2+10*c+2]+I*y[2+10*c+3];
		xi20=y[2+10*c+4]+I*y[2+10*c+5];
		xi21=y[2+10*c+6]+I*y[2+10*c+7];
		xi22=y[2+10*c+8]+I*y[2+10*c+9];

		dxi2m2dt=-I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi2m1;
		dxi2m1dt=-I*sqrt(6*L*(L+1))*timephase(2.0f,t,config)*xi20-I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi2m2;

 		if(c!=0)
		{
			double complex invphase2m1;

			invphase2m1=timephase(-(L*(L+1.0f)+omegak(k)-2.0f),t,config);

			dxi2m1dt+=I*6.0f*V2(k,localdensity,config)/omegak(k)*invphase2m1*localAminus;
			dxi2m1dt+=I*V2(k,localdensity,config)/omegak(k)*sqrt(6*L*(L+1))*invphase2m1*g;
		}

		dxi20dt=-I*sqrt(6*L*(L+1))*timephase(-2.0f,t,config)*(xi2m1+xi21);
		dxi21dt=-I*sqrt(6*L*(L+1))*timephase(2.0f,t,config)*xi20-I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi22;

 		if(c!=0)
		{
			double complex invphase21;

			invphase21=timephase(-(L*(L+1.0f)+omegak(k)-2.0f),t,config);

			dxi21dt+=I*6.0f*V2(k,localdensity,config)/omegak(k)*invphase21*localAplus;
			dxi21dt+=I*V2(k,localdensity,config)/omegak(k)*sqrt(6*L*(L+1))*invphase21*g;
		}

		dxi22dt=-I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi21;

		dydt[2+10*c]=creal(dxi2m2dt);
		dydt[2+10*c+1]=cimag(dxi2m2dt);
		dydt[2+10*c+2]=creal(dxi2m1dt);
		dydt[2+10*c+3]=cimag(dxi2m1dt);
		dydt[2+10*c+4]=creal(dxi20dt);
		dydt[2+10*c+5]=cimag(dxi20dt);
		dydt[2+10*c+6]=creal(dxi21dt);
		dydt[2+10*c+7]=cimag(dxi21dt);
		dydt[2+10*c+8]=creal(dxi22dt);
		dydt[2+10*c+9]=cimag(dxi22dt);
	}

	return GSL_SUCCESS;
}
