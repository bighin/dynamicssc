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
#include "coherent.h"

#define MIN(x,y)		(((x)<(y))?(x):(y))

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

double omegak(double k,struct configuration_t *config)
{
	double num,den;
	double a0,a1,a2,a3,a4,a5,b0,b1,b2,b3,b4,b5;

	if(config->dispersion_is_experimental==false)
		return k*config->soundspeed;

	/*
		Padè approximant for the dispersion relation
	*/

	if(config->moleculetype==MOLECULE_I2)
	{
		/*
			Coefficients for I2
		*/

		a0=0.0f;
		a1=78979.82692539008;
		a2=-4453.106390752133;
		a3=383.5583835895215;
		a4=-21.92511621355489;
		a5=0.3824072227009381;

		b0=3433.7252161679503;
		b1=-218.02760941323945;
		b2=12.467652732020886;
		b3=0.5009420138478766;
		b4=-0.05760613914422612;
		b5=0.0010925920648598231;
	}
	else if(config->moleculetype==MOLECULE_CS2)
	{
		/*
			Coefficients for CS2
		*/

		a0=0.0f;
		a1=39815.65200133856;
		a2=-4031.0062614672765;
		a3=570.0463741407423;
		a4=-54.18214205687705;
		a5=1.6053060335262601;

		b0=2961.928011862246;
		b1=-334.1600104773988;
		b2=32.52943513283185;
		b3=1.9844702294478243;
		b4=-0.409499213901235;
		b5=0.013377550279385501;
	}
	else if(config->moleculetype==MOLECULE_OCS)
	{
		/*
			Coefficients for OCS
		*/

		a0=0.0f;
		a1=22655.7;
		a2=-4009.93;
		a3=602.768;
		a4=-62.1808;
		a5=2.3052;

		b0=2256.28;
		b1=-421.22;
		b2=53.0869;
		b3=1.39869;
		b4=-0.713085;
		b5=0.0329314;
	}
	else
	{
		fprintf(stderr,"Unknown molecule type!\n");
		exit(0);
	}

	num=a0+k*(a1+k*(a2+k*(a3+k*(a4+k*a5))));
	den=b0+k*(b1+k*(b2+k*(b3+k*(b4+k*b5))));

	return num/den;
}

/*
	Potentials
*/

double U0(double n,double k,struct configuration_t *config)
{
	double a,b,c;

	a=sqrt((8.0f*n*k*k*ek(k))/(omegak(k,config)));
	b=exp(-0.5f*k*k*pow(config->r0,2.0));
	c=4*M_PI*pow(config->r0,-3.0f);

	return config->u0*a*b/c;
}

double U2gaussian(double n,double k,struct configuration_t *config)
{
	double a,b,c,d,z;

	z=k*config->r2;
	a=sqrt((8.0f*n*k*k*ek(k))/(5.0f*omegak(k,config)));
	b=-2.0f*exp(-0.5*z*z)*z*(3.0f+z*z);
	c=3.0f*sqrt(2.0f*M_PI)*erf(z/sqrt(2.0f));
	d=8.0f*pow(k,3.0f)*M_PI;
	
	return config->u2*a*(b+c)/d;
}

double U2morse(double n,double k,struct configuration_t *config)
{
	double A=1.25f;
	double re=0.75f;

	double x1,x2,x3,x4,x5,x6,z;

	double A2=pow(A,2.0f);
	double A3=pow(A,3.0f);
	double k2=pow(k,2.0f);
	double k3=pow(k,3.0f);

	double expAre=exp(A*re);

	x1=(6*A3*k)/pow(A2+k2,2.0f);
	x2=(10.0f*A*k3)/pow(A2+k2,2.0f);
	x3=-(24.0f*A3*expAre*k)/pow(4.0f*A2+k2,2.0f);
	x4=-(10.0f*A*expAre*k3)/pow(4.0f*A2+k2,2.0f);

	x5=3.0f*expAre*arccot(2.0f*A/k);
	x6=-6.0f*atan(k/A);

	z=-expAre/(2.0f*sqrt(2.0f)*k3*pow(M_PI,3.0f/2.0f));
	
	return config->u2*sqrt((8.0f*n*k*k*ek(k))/(5.0f*omegak(k,config)))*z*(x1+x2+x3+x4+x5+x6);
}

double U2(double n,double k,struct configuration_t *config)
{
	if(config->morse==true)
		return U2morse(n,k,config);

	return U2gaussian(n,k,config);
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

double complex Pk(double complex En,int L,double k,struct configuration_t *config)
{
        double complex a,b,c;

        a=-En+L*(L+1)+omegak(k,config)+4.0f;
        b=-En+L*(L+1)+omegak(k,config)+6.0f;
        c=-En+L*(L+1)+omegak(k,config)-2.0f;

	return a-(12.0f*L*(L+1))/b-(4.0f*(L*(L+1)-2.0f))/c;
}

double complex fscale(double k,int L,struct configuration_t *config)
{
	double epsilon=0.001f;
	double en=L*(L+1);

	if(config->fscale==true)
		return Pk(0.9*en+I*epsilon,L,k,config)/(1.0f+k*k)+I*epsilon;

	return 1.0;
}

int fnorm(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the 'external' variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables */

        double k;
	double complex xi2m2,xi2m1,xi20,xi21,xi22;
	int L;

	k=x[0];

	xi2m2=get_point(container->intrexi2m2,k)+I*get_point(container->intimxi2m2,k);
	xi2m1=get_point(container->intrexi2m1,k)+I*get_point(container->intimxi2m1,k);
	xi20=get_point(container->intrexi20,k)+I*get_point(container->intimxi20,k);
	xi21=get_point(container->intrexi21,k)+I*get_point(container->intimxi21,k);
	xi22=get_point(container->intrexi22,k)+I*get_point(container->intimxi22,k);

	L=container->params->L;

	xi2m2/=fscale(k,L,config);
	xi2m1/=fscale(k,L,config);
	xi20/=fscale(k,L,config);
	xi21/=fscale(k,L,config);
	xi22/=fscale(k,L,config);

	fval[0]=conj(xi2m2)*xi2m2+conj(xi2m1)*xi2m1+conj(xi20)*xi20+conj(xi21)*xi21+conj(xi22)*xi22;

	return 0;
}

double norm_phonons_1phonon(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	double *x,*yre2m2,*yim2m2,*yre2m1,*yim2m1,*yre20,*yim20,*yre21,*yim21,*yre22,*yim22;
	double xmin,xmax,err;

	struct container_t container;

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
                yre2m2[c]=y[10+10*c];
                yim2m2[c]=y[10+10*c+1];
                yre2m1[c]=y[10+10*c+2];
                yim2m1[c]=y[10+10*c+3];
                yre20[c]=y[10+10*c+4];
                yim20[c]=y[10+10*c+5];
                yre21[c]=y[10+10*c+6];
                yim21[c]=y[10+10*c+7];
                yre22[c]=y[10+10*c+8];
                yim22[c]=y[10+10*c+9];
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

	container.params=params;

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

double norm_phonons_free(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	return 0.0f;
}

double norm_phonons(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	switch(config->evolution)
	{
		case EVOLUTION_FREE:
		return norm_phonons_free(t,y,params,config);

		case EVOLUTION_1PHONON:
		case EVOLUTION_1PHONONFT:
		return norm_phonons_1phonon(t,y,params,config);

		case EVOLUTION_COHERENT:
		return norm_phonons_coherent(t,y,params,config);
	}

	fprintf(stderr,"Unkown evolution type!\n");
	exit(0);

	/*
		Never reached.
	*/

	return 0.0f;
}

double norm_qp_free(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	double complex g=y[0]+I*y[1];

	/*
		Here we should multiply g by a phase, but since we are just interested
		in the modulus, this is not needed...
	*/
	
	return sqrt(conj(g)*g);
}

double norm_qp_1phonon(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	double complex g=y[0]+I*y[1];

	/*
		Here we should multiply g by a phase, but since we are just interested
		in the modulus, this is not needed...
	*/
	
	return sqrt(conj(g)*g);
}

double norm_qp(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	switch(config->evolution)
	{
		case EVOLUTION_FREE:
		return norm_qp_free(t,y,params,config);

		case EVOLUTION_1PHONON:
		case EVOLUTION_1PHONONFT:
		return norm_qp_1phonon(t,y,params,config);

		case EVOLUTION_COHERENT:
		return norm_qp_coherent(t,y,params,config);
	}

	fprintf(stderr,"Unkown evolution type!\n");
	exit(0);

	/*
		Never reached.
	*/

	return 0.0f;
}

double norm(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	switch(config->evolution)
	{
		case EVOLUTION_FREE:
		case EVOLUTION_1PHONON:
		case EVOLUTION_1PHONONFT:
		return sqrt(pow(norm_qp(t,y,params,config),2.0f)+
		            pow(norm_phonons(t,y,params,config),2.0f));

		case EVOLUTION_COHERENT:
		return norm_coherent(t,y,params,config);
	}

	fprintf(stderr,"Unkown evolution type!\n");
	exit(0);

	/*
		Never reached.
	*/

	return 0.0f;
}

double W(double k,struct configuration_t *config)
{
	if(config->wtype==1)
		return omegak(k,config);

	if(config->wtype==2)
		return omegak(k,config)+6.0f;

	fprintf(stderr,"Fatal error: wrong wtype in configuration");
	exit(0);

	return 0.0f;
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

	phase21=timephase(-(omegak(k,config)+4.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase21*(get_point(container->intrexi21,k)+I*get_point(container->intimxi21,k));
	f/=fscale(k,L,config);

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
                yre[c]=y[10+10*c+6];
                yim[c]=y[10+10*c+7];
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

	phase2m1=timephase(-(omegak(k,config)+4.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase2m1*(get_point(container->intrexi2m1,k)+I*get_point(container->intimxi2m1,k));
	f/=fscale(k,L,config);

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
                yre[c]=y[10+10*c+2];
                yim[c]=y[10+10*c+3];
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

int fB(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase20;
	double t,localdensity;
	int L;

	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;

	k=x[0];

	phase20=timephase(-(omegak(k,config)+6.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase20*(get_point(container->intrexi20,k)+I*get_point(container->intimxi20,k))*(omegak(k,config)+6.0-W(k,config));
	f/=fscale(k,L,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex B(double t,const double y[],struct params_t *params,double localdensity)
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
                yre[c]=y[10+10*c+4];
                yim[c]=y[10+10*c+5];
	}

        container.intrexi20=init_interpolation(x,yre,config->gridpoints);
        container.intimxi20=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fB,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi20);
	fini_interpolation(container.intimxi20);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int sc_time_evolution_1phonon(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	double complex localAplus,localAminus,localB;
	double complex g,dgdt;
	double localdensity;
	int c,L;

	L=params->L;

	g=y[0]+I*y[1];
	dgdt=0.0f;

        if(config->centrifugal==true)
                dgdt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*g;

	if(config->ramp==true)
		localdensity=config->density*adiabatic_ramp(t,config);
	else
		localdensity=config->density;

	localAplus=Aplus(t,y,p,localdensity);
	localAminus=Aminus(t,y,p,localdensity);
	localB=B(t,y,p,localdensity);

	dgdt+=-I*sqrt(6*L*(L+1))*(localAplus+localAminus);
	dgdt+=I*localB;

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

		xi2m2=y[10+10*c]+I*y[10+10*c+1];
		xi2m1=y[10+10*c+2]+I*y[10+10*c+3];
		xi20=y[10+10*c+4]+I*y[10+10*c+5];
		xi21=y[10+10*c+6]+I*y[10+10*c+7];
		xi22=y[10+10*c+8]+I*y[10+10*c+9];

		dxi2m2dt=I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi2m1;
		dxi2m1dt=I*sqrt(6*L*(L+1))*timephase(-2.0f,t,config)*xi20+I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi2m2;

 		if(c!=0)
		{
			double complex invphase2m1;

			invphase2m1=timephase(omegak(k,config)+4.0f,t,config);

			dxi2m1dt+=-I*6.0f*V2(k,localdensity,config)/W(k,config)*invphase2m1*fscale(k,L,config)*localAminus;
			dxi2m1dt+=-I*V2(k,localdensity,config)/W(k,config)*sqrt(6*L*(L+1))*invphase2m1*fscale(k,L,config)*g;
		}

		dxi20dt=I*sqrt(6*L*(L+1))*timephase(2.0f,t,config)*(xi2m1+xi21);

		if(c!=0)
			dxi20dt+=I*g*V2(k,localdensity,config)/W(k,config)*(omegak(k,config)+6.0-W(k,config))*timephase(omegak(k,config)+6.0f,t,config)*fscale(k,L,config);

		dxi21dt=I*sqrt(6*L*(L+1))*timephase(-2.0f,t,config)*xi20+I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi22;

 		if(c!=0)
		{
			double complex invphase21;

			invphase21=timephase(omegak(k,config)+4.0f,t,config);

			dxi21dt+=-I*6.0f*V2(k,localdensity,config)/W(k,config)*invphase21*fscale(k,L,config)*localAplus;
			dxi21dt+=-I*V2(k,localdensity,config)/W(k,config)*sqrt(6*L*(L+1))*invphase21*fscale(k,L,config)*g;
		}

		dxi22dt=I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi21;

                if(config->centrifugal==true)
                {
			dxi2m2dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi2m2;
			dxi2m1dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi2m1;
			dxi20dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi20;
			dxi21dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi21;
			dxi22dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi22;
		}

		dydt[10+10*c]=creal(dxi2m2dt);
		dydt[10+10*c+1]=cimag(dxi2m2dt);
		dydt[10+10*c+2]=creal(dxi2m1dt);
		dydt[10+10*c+3]=cimag(dxi2m1dt);
		dydt[10+10*c+4]=creal(dxi20dt);
		dydt[10+10*c+5]=cimag(dxi20dt);
		dydt[10+10*c+6]=creal(dxi21dt);
		dydt[10+10*c+7]=cimag(dxi21dt);
		dydt[10+10*c+8]=creal(dxi22dt);
		dydt[10+10*c+9]=cimag(dxi22dt);
	}

	return GSL_SUCCESS;
}

int sc_time_evolution_free(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	double complex g,dgdt;
	int L=params->L;

	g=y[0]+I*y[1];
	dgdt=0.0f;

        if(config->centrifugal==true)
                dgdt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*g;

	dydt[0]=creal(dgdt);
	dydt[1]=cimag(dgdt);

	return GSL_SUCCESS;
}

/*
	The Bose-Einstein distribution
*/

double fbe(struct configuration_t *config,double x)
{
	double temperature,beta,kelvins;

	kelvins=config->temperature;

	if(config->moleculetype==MOLECULE_I2)
	{
		/*
			Inverse temperature of the bath, in units of B, for I2
		*/

		temperature=kelvins/0.0536459;
		beta=1.0f/temperature;
	}
	else if(config->moleculetype==MOLECULE_CS2)
	{
		/*
			Inverse temperature of the bath, in units of B, for CS2
		*/

		temperature=kelvins/0.156839;
		beta=1.0f/temperature;
	}
	else if(config->moleculetype==MOLECULE_OCS)
	{
		/*
			Inverse temperature of the bath, in units of B, for OCS
		*/

		temperature=kelvins/0.291866;
		beta=1.0f/temperature;
	}
	else
	{
		fprintf(stderr,"Unknown molecule type!\n");
		exit(0);
	}

	return 1.0f/(exp(beta*x)-1.0f);
}

int fAplus_thermal(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase21;
	double t,localdensity,fbek;
	int L;
	
	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;

	k=x[0];
	
	/*
		Note that here we can take fbek to be zero at zero,
		since the macroscopic occupation of the ground state
		has been accounted for a long time ago while doing
		Bogoliubov transformation.

		The same applies to fAplus_thermal() and to fB_thermal().

		Probably this check is not even needed, since
		the integration routine shouldn't evaluate the integral
		at the edges of the integration domain.
	*/

#define ALMOST_ZERO(x)	(fabs(x)<=1e-8)

	fbek=(ALMOST_ZERO(k))?(0.0f):(fbe(config,omegak(k,config)));

	phase21=timephase(-(omegak(k,config)+4.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase21*(get_point(container->intrexi21,k)+I*get_point(container->intimxi21,k));
	f*=(1.0f+fbek);
	f/=fscale(k,L,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex Aplus_thermal(double t,const double y[],struct params_t *params,double localdensity)
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
                yre[c]=y[10+10*c+6];
                yim[c]=y[10+10*c+7];
	}

        container.intrexi21=init_interpolation(x,yre,config->gridpoints);
        container.intimxi21=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fAplus_thermal,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi21);
	fini_interpolation(container.intimxi21);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int fAminus_thermal(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase2m1;
	double t,localdensity,fbek;
	int L;
	
	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;

	k=x[0];
	fbek=(ALMOST_ZERO(k))?(0.0f):(fbe(config,omegak(k,config)));

	phase2m1=timephase(-(omegak(k,config)+4.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase2m1*(get_point(container->intrexi2m1,k)+I*get_point(container->intimxi2m1,k));
	f*=(1.0f+fbek);
	f/=fscale(k,L,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex Aminus_thermal(double t,const double y[],struct params_t *params,double localdensity)
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
                yre[c]=y[10+10*c+2];
                yim[c]=y[10+10*c+3];
	}

        container.intrexi2m1=init_interpolation(x,yre,config->gridpoints);
        container.intimxi2m1=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fAminus_thermal,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi2m1);
	fini_interpolation(container.intimxi2m1);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int fB_thermal(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase20;
	double t,localdensity,fbek;
	int L;

	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;

	k=x[0];
	fbek=(ALMOST_ZERO(k))?(0.0f):(fbe(config,omegak(k,config)));

	phase20=timephase(-(omegak(k,config)+6.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase20*(get_point(container->intrexi20,k)+I*get_point(container->intimxi20,k))*(omegak(k,config)+6.0-W(k,config));
	f*=(1.0f+fbek);
	f/=fscale(k,L,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex B_thermal(double t,const double y[],struct params_t *params,double localdensity)
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
                yre[c]=y[10+10*c+4];
                yim[c]=y[10+10*c+5];
	}

        container.intrexi20=init_interpolation(x,yre,config->gridpoints);
        container.intimxi20=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fB_thermal,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi20);
	fini_interpolation(container.intimxi20);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int fC_thermal(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct container_t *container=(struct container_t *)(fdata);
	struct configuration_t *config=container->params->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f,phase20;
	double t,localdensity,fbek;
	int L;

	t=container->t;
	localdensity=container->localdensity;
	L=container->params->L;

	k=x[0];
	fbek=(ALMOST_ZERO(k))?(0.0f):(fbe(config,omegak(k,config)));

	phase20=timephase(-(omegak(k,config)+6.0f),t,config);
	f=V2(k,localdensity,config)/W(k,config)*phase20*(get_point(container->intrexi20,k)+I*get_point(container->intimxi20,k));
	f*=(1.0f+fbek);
	f*=fbek;
	f/=fscale(k,L,config);

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

double complex C_thermal(double t,const double y[],struct params_t *params,double localdensity)
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
                yre[c]=y[10+10*c+4];
                yim[c]=y[10+10*c+5];
	}

        container.intrexi20=init_interpolation(x,yre,config->gridpoints);
        container.intimxi20=init_interpolation(x,yim,config->gridpoints);
	container.t=t;
	container.params=params;
	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fC_thermal,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intrexi20);
	fini_interpolation(container.intimxi20);

	if(x)	free(x);
	if(yre)	free(yre);
	if(yim)	free(yim);

	return res[0]+I*res[1];
}

int sc_time_evolution_1phononft(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	double complex localAplus,localAminus,localB,localC;
	double complex g,dgdt;
	double localdensity;
	int c,L;

	L=params->L;

	g=y[0]+I*y[1];

	dgdt=0.0f;

        if(config->centrifugal==true)
                dgdt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*g;

	if(config->ramp==true)
		localdensity=config->density*adiabatic_ramp(t,config);
	else
		localdensity=config->density;

	localAplus=Aplus_thermal(t,y,p,localdensity);
	localAminus=Aminus_thermal(t,y,p,localdensity);
	localB=B_thermal(t,y,p,localdensity);
	localC=C_thermal(t,y,p,localdensity);

	dgdt+=-I*sqrt(6*L*(L+1))*(localAplus+localAminus);
	dgdt+=I*localB;

#warning Check the next line! (g phase and alpha phase)
	dgdt+=-I*12.0f*localC;

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

		/*
			The Bose-Einstein distribution for omega(k).
		
			Note that is does not make sense for c=0, i.e. k=0.
		*/

		double fbek=(c==0)?(0.0f):(fbe(config,omegak(k,config)));

		xi2m2=y[10+10*c]+I*y[10+10*c+1];
		xi2m1=y[10+10*c+2]+I*y[10+10*c+3];
		xi20=y[10+10*c+4]+I*y[10+10*c+5];
		xi21=y[10+10*c+6]+I*y[10+10*c+7];
		xi22=y[10+10*c+8]+I*y[10+10*c+9];

		dxi2m2dt=I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi2m1*(1.0f+fbek);
		dxi2m1dt=I*sqrt(6*L*(L+1))*timephase(-2.0f,t,config)*xi20*(1.0f+fbek)+I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi2m2*(1.0f+fbek);

 		if(c!=0)
		{
			double complex invphase2m1;

			invphase2m1=timephase(omegak(k,config)+4.0f,t,config);

			dxi2m1dt+=-I*6.0f*V2(k,localdensity,config)/W(k,config)*invphase2m1*fscale(k,L,config)*localAminus;
			dxi2m1dt+=-I*V2(k,localdensity,config)/W(k,config)*sqrt(6*L*(L+1))*invphase2m1*fscale(k,L,config)*g;
		}

		dxi20dt=I*sqrt(6*L*(L+1))*timephase(2.0f,t,config)*(xi2m1+xi21)*(1.0f+fbek);

		if(c!=0)
		{
			dxi20dt+=I*g*V2(k,localdensity,config)/W(k,config)*(omegak(k,config)+6.0-W(k,config))*timephase(omegak(k,config)+6.0f,t,config)*fscale(k,L,config);

#warning Check the next line (g phase and alpha phase)
		        dxi20dt+=-12.0f*I*g*V2(k,localdensity,config)/W(k,config)*fbek*timephase(omegak(k,config)+6.0f,t,config)*fscale(k,L,config);
		}

		dxi21dt=I*sqrt(6*L*(L+1))*timephase(-2.0f,t,config)*xi20*(1.0f+fbek)+I*2.0f*sqrt(L*(L+1)-2)*timephase(6.0f,t,config)*xi22*(1.0f+fbek);

 		if(c!=0)
		{
			double complex invphase21;

			invphase21=timephase(omegak(k,config)+4.0f,t,config);

			dxi21dt+=-I*6.0f*V2(k,localdensity,config)/W(k,config)*invphase21*fscale(k,L,config)*localAplus;
			dxi21dt+=-I*V2(k,localdensity,config)/W(k,config)*sqrt(6*L*(L+1))*invphase21*fscale(k,L,config)*g;
		}

		dxi22dt=I*2.0f*sqrt(L*(L+1)-2)*timephase(-6.0f,t,config)*xi21*(1.0f+fbek);

		/*
			Additional phases which are not completely removed by the transformation.

			We could include them in the transformation, but that would mean using two
			different transformations at T=0 and at finite temperature.
		
			For simplicity's sake, since these phases are slowly-oscillating, we leave
			them here.
		*/

		dxi2m2dt+=-I*(omegak(k,config)-2.0f)*fbek*xi2m2;
		dxi2m1dt+=-I*(omegak(k,config)+4.0f)*fbek*xi2m1;
		dxi20dt+=-I*(omegak(k,config)+6.0f)*fbek*xi20;
		dxi21dt+=-I*(omegak(k,config)+4.0f)*fbek*xi21;
		dxi22dt+=-I*(omegak(k,config)-2.0f)*fbek*xi22;

                if(config->centrifugal==true)
                {
			dxi2m2dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi2m2;
			dxi2m1dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi2m1;
			dxi20dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi20;
			dxi21dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi21;
			dxi22dt+=I*config->centrifugalD*pow(MIN(L,config->centrifugalLcutoff)*(MIN(L,config->centrifugalLcutoff)+1),2.0f)*xi22;
		}

		dydt[10+10*c]=creal(dxi2m2dt);
		dydt[10+10*c+1]=cimag(dxi2m2dt);
		dydt[10+10*c+2]=creal(dxi2m1dt);
		dydt[10+10*c+3]=cimag(dxi2m1dt);
		dydt[10+10*c+4]=creal(dxi20dt);
		dydt[10+10*c+5]=cimag(dxi20dt);
		dydt[10+10*c+6]=creal(dxi21dt);
		dydt[10+10*c+7]=cimag(dxi21dt);
		dydt[10+10*c+8]=creal(dxi22dt);
		dydt[10+10*c+9]=cimag(dxi22dt);
	}

	return GSL_SUCCESS;
}

int sc_time_evolution(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	switch(config->evolution)
	{
		case EVOLUTION_FREE:
		return sc_time_evolution_free(t,y,dydt,p);

		case EVOLUTION_1PHONON:
		return sc_time_evolution_1phonon(t,y,dydt,p);

		case EVOLUTION_1PHONONFT:
		return sc_time_evolution_1phononft(t,y,dydt,p);

		case EVOLUTION_COHERENT:
		return sc_time_evolution_coherent(t,y,dydt,p);
	
		default:
		fprintf(stderr,"Unkown evolution type!\n");
		exit(0);
	}

	/*
		Never reached
	*/

	return GSL_FAILURE;
}
