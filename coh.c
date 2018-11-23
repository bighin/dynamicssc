#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>

#include "coh.h"
#include "dsc.h"
#include "bigpsi.h"
#include "observables.h"

double complex Gminus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret;
	int offsetL;

	ret=0.0f;
	offsetL=L*(10+10*config->gridpoints);

	for(int n=-2;n<=2;n++)
	{
		int nminus;
		double complex gn,gnminus;
		double Cminus;

		nminus=n-1;

		if(abs(nminus)>2)
			continue;

		gn=psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1];
		gnminus=psi->y[offsetL+2*(nminus+2)]+I*psi->y[offsetL+2*(nminus+2)+1];
		Cminus=sqrtf(L*(L+1.0f))*cg(L,nminus,1,1,L,n);

		ret+=conj(gn)*gnminus*Cminus;
	}

	return ret;
}

double complex Gplus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret;
	int offsetL;
	
	ret=0.0f;
	offsetL=L*(10+10*config->gridpoints);

	for(int n=-2;n<=2;n++)
	{
		int nplus;
		double complex gn,gnplus;
		double Cplus;

		nplus=n+1;

		if(abs(nplus)>2)
			continue;

		gn=psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1];
		gnplus=psi->y[offsetL+2*(nplus+2)]+I*psi->y[offsetL+2*(nplus+2)+1];
		Cplus=sqrtf(L*(L+1.0f))*cg(L,nplus,1,-1,L,n);

		ret+=conj(gn)*gnplus*Cplus;
	}

	return ret;
}

double complex Gzero(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret;
	int offsetL;
	
	ret=0.0f;
	offsetL=L*(10+10*config->gridpoints);

	for(int n=-2;n<=2;n++)
	{
		double complex gn;

		gn=psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1];

		ret+=n*conj(gn)*gn;
	}

	return ret;
}

double complex Lambdaplus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int plus=1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(plus,mu,nu)*Dcross(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdaminus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int minus=1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(minus,mu,nu)*Dcross(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdazero(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int zero=1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(zero,mu,nu)*Dcross(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Vbeta(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	return creal(Dsingle(psi,L,L,0,DINT_MODE_SUPERPLAIN_VK0,config));
}

int sc_time_evolution_coh(double t,const double y[],double dydt[],void *p)
{
	struct params_t *params=(struct params_t *)(p);
	struct configuration_t *config=params->config;

	struct bigpsi_t *psi=(struct bigpsi_t *)(params->parent);

	double complex lambdaplus,lambdaminus,lambdazero;
	double complex gplus,gminus,gzero;

	int plus,minus,zero;

	int L=params->L;

	plus=+1;
	minus=-1;
	zero=0;

	gplus=Gplus(psi,L,config);
	gminus=Gminus(psi,L,config);
	gzero=Gzero(psi,L,config);

	lambdaplus=Lambdaplus(psi,L,config);
	lambdaminus=Lambdaminus(psi,L,config);
	lambdazero=Lambdazero(psi,L,config);

	if((config->ramp==true)||(config->bosonsfinitetemperature==true)||(config->centrifugal))
	{
		printf("Some features (adiabatic ramp, finite temperature evolution, centrifugal distortion) are not supported at the moment when using coherent state evolution.\n");
		exit(0);
	}

	if(config->freeevolution==true)
	{
		for(int n=-2;n<=2;n++)
		{
			dydt[2*(n+2)]=0.0f;
			dydt[2*(n+2)+1]=0.0f;
		}

		return GSL_SUCCESS;
	}

	for(int n=-2;n<=2;n++)
	{
		double complex gn,dgndt;
		double complex lambdasquared;
		int nplus,nminus;

		gn=y[2*(n+2)]+I*y[2*(n+2)+1];
		dgndt=0.0f;

		lambdasquared=conj(lambdaplus)*lambdaplus+
			      conj(lambdaminus)*lambdaminus+
			      conj(lambdazero)*lambdazero;

		dgndt+=-I*gn*Vbeta(psi,L,config);
		dgndt+=I*gn*lambdasquared;

		dgndt+=-I*gzero*lambdazero;
		dgndt+=-I*gplus*lambdaminus;
		dgndt+=-I*gminus*lambdaplus;

		dgndt+=2.0f*I*n*gn*lambdazero;

		nminus=n-1;
		if(abs(nminus)<=2)
		{
			double complex gnminus=y[2*(nminus+2)]+I*y[2*(nminus+2)+1];
			double Cminus=sqrtf(L*(L+1.0f))*cg(L,nminus,1,1,L,n);

			dgndt+=2.0f*I*gnminus*Cminus*lambdaplus;
		}

		nplus=n+1;
		if(abs(nplus)<=2)
		{
			double complex gnplus=y[2*(nplus+2)]+I*y[2*(nplus+2)+1];
			double Cplus=sqrtf(L*(L+1.0f))*cg(L,nplus,1,-1,L,n);

			dgndt+=2.0f*I*gnplus*Cplus*lambdaminus;
		}
		
		dydt[2*(n+2)]=creal(dgndt);
		dydt[2*(n+2)+1]=cimag(dgndt);
	}

	if(L<=0)
		return GSL_SUCCESS;

	for(int c=0;c<config->gridpoints;c++)
	{
		double gridstep=config->cutoff/config->gridpoints;
		double k=gridstep*c;

		double complex xi2[5];

		for(int mu=-2;mu<=2;mu++)
			xi2[2+mu]=y[10+10*c+(mu+2)*2]+I*y[10+10*c+(mu+2)*2+1];

		for(int mu=-2;mu<=2;mu++)
		{
			double complex dxi2mudt=0.0f;

			if(mu==0)
				dxi2mudt+=-I*V2(k,config->density,config);

#warning Fattori di fase sui vari beta!!!

			dxi2mudt+=-I*(omegak(k,config)+6.0-2.0f*mu*gzero)*xi2[2+mu];

			for(int nu=-2;nu<=2;nu++)
			{
				double complex prod=conj(lambdaplus)*sigma_matrix(plus,mu,nu)+
				                    conj(lambdaminus)*sigma_matrix(minus,mu,nu)+
				                    conj(lambdazero)*sigma_matrix(zero,mu,nu);

				dxi2mudt+=-2.0f*I*prod*xi2[2+nu];
			}

			for(int nu=-2;nu<=2;nu++)
			{
				double complex prod=gplus*sigma_matrix(minus,mu,nu)+
				                    gminus*sigma_matrix(plus,mu,nu);

				dxi2mudt+=2.0f*I*prod*xi2[2+nu];
			}

			dydt[10+10*c+8]=creal(dxi2mudt);
			dydt[10+10*c+9]=cimag(dxi2mudt);	
		}
	}

	return GSL_SUCCESS;
}

int big_sc_time_evolution_coh(double t,const double y[],double dydt[],void *data)
{
	struct bigpsi_t *bigpsi=data;
	struct params_t *params=bigpsi->params;
	struct configuration_t *config=params->config;

	double complex C;

	for(int c=0;c<(10+10*config->gridpoints)*bigpsi->nrpsis;c++)
		dydt[c]=0.0f;

	/*
		We apply the strong coupling evolution to each different L state...
	*/

	for(int c=0;c<bigpsi->nrpsis;c++)
	{
		int offset=c*(10+10*config->gridpoints);

		sc_time_evolution_coh(t,&y[offset],&dydt[offset],&params[c]);
	}

#if 0
	/*
		and then the level mixing due to the electric field!
	*/

	if(laser_is_on(config,t)==false)
		return GSL_SUCCESS;

	C=(2.0f/3.0f)*get_laser_intensity(config->fluence,config->duration,t,config);

#define GETL(x) (params[x].L)
#define GETM() 	(params[0].M)

	for(int d=0;d<bigpsi->nrpsis;d++)
	{
		int offset=d*(10+10*config->gridpoints);
		double complex gLM,dgLMdt;
		int L=GETL(d);

		gLM=(y[offset+0]+I*y[offset+1]);

		dgLMdt=0.0f;
		for(int c=d-2;c<=d+2;c++)
		{
			int local_offset=c*(10+10*config->gridpoints);
			double complex gLprimeM;
			int Lprime;

			if((c<0)||(c>=bigpsi->nrpsis))
				continue;

			Lprime=GETL(c);
			gLprimeM=y[local_offset+0]+I*y[local_offset+1];

			dgLMdt+=I*C*Q(GETL(d),GETL(c),GETM(),0)*gLprimeM*timephase((L*(L+1.0f)-Lprime*(Lprime+1.0f)),t,config);
		}

		dgLMdt+=I*(C/2.0f)*gLM;

		dydt[offset+0]+=creal(dgLMdt);
		dydt[offset+1]+=cimag(dgLMdt);

		/*
			The phonons with L=0 are not coupled to anything!
		*/

		if(L==0)
			continue;

		for(int i=0;i<config->gridpoints;i++)
		{
			double complex xiLM2m2,xiLM2m1,xiLM20,xiLM21,xiLM22;
			double complex dxiLM2m2dt,dxiLM2m1dt,dxiLM20dt,dxiLM21dt,dxiLM22dt;
			int L=GETL(d);
			
			double gridstep=config->cutoff/config->gridpoints;
			double k=i*gridstep;

			dxiLM2m2dt=dxiLM2m1dt=dxiLM20dt=dxiLM21dt=dxiLM22dt=0.0f;
			for(int c=d-2;c<=d+2;c++)
			{
				int local_offset=c*(10+10*config->gridpoints);

				double complex xiLprimeM2m2,xiLprimeM2m1,xiLprimeM20,xiLprimeM21,xiLprimeM22;
				double complex phasediff;
				double complex scalefactor;

				int Lprime;

				if((c<1)||(c>=bigpsi->nrpsis))
					continue;

				Lprime=GETL(c);

				phasediff=timephase((L*(L+1.0f)-Lprime*(Lprime+1.0f)),t,config);

				xiLprimeM2m2=y[local_offset+10+10*i]+I*y[local_offset+10+10*i+1];
				xiLprimeM2m1=y[local_offset+10+10*i+2]+I*y[local_offset+10+10*i+3];
				xiLprimeM20=y[local_offset+10+10*i+4]+I*y[local_offset+10+10*i+5];
				xiLprimeM21=y[local_offset+10+10*i+6]+I*y[local_offset+10+10*i+7];
				xiLprimeM22=y[local_offset+10+10*i+8]+I*y[local_offset+10+10*i+9];

				scalefactor=fscale(k,L,config)/fscale(k,Lprime,config);

				dxiLM2m2dt+=I*phasediff*C*Q(L,Lprime,GETM(),-2)*xiLprimeM2m2*scalefactor;
				dxiLM2m1dt+=I*phasediff*C*Q(L,Lprime,GETM(),-1)*xiLprimeM2m1*scalefactor;
				dxiLM20dt+=I*phasediff*C*Q(L,Lprime,GETM(),0)*xiLprimeM20*scalefactor;
				dxiLM21dt+=I*phasediff*C*Q(L,Lprime,GETM(),1)*xiLprimeM21*scalefactor;
				dxiLM22dt+=I*phasediff*C*Q(L,Lprime,GETM(),2)*xiLprimeM22*scalefactor;
			}

			xiLM2m2=y[offset+10+10*i]+I*y[offset+10+10*i+1];
			xiLM2m1=y[offset+10+10*i+2]+I*y[offset+10+10*i+3];
			xiLM20=y[offset+10+10*i+4]+I*y[offset+10+10*i+5];
			xiLM21=y[offset+10+10*i+6]+I*y[offset+10+10*i+7];
			xiLM22=y[offset+10+10*i+8]+I*y[offset+10+10*i+9];

			dxiLM2m2dt+=I*(C/2.0f)*xiLM2m2;
			dxiLM2m1dt+=I*(C/2.0f)*xiLM2m1;
			dxiLM20dt+=I*(C/2.0f)*xiLM20;
			dxiLM21dt+=I*(C/2.0f)*xiLM21;
			dxiLM22dt+=I*(C/2.0f)*xiLM22;

			dydt[offset+10+10*i]+=creal(dxiLM2m2dt);
			dydt[offset+10+10*i+1]+=cimag(dxiLM2m2dt);

			dydt[offset+10+10*i+2]+=creal(dxiLM2m1dt);
			dydt[offset+10+10*i+3]+=cimag(dxiLM2m1dt);

			dydt[offset+10+10*i+4]+=creal(dxiLM20dt);
			dydt[offset+10+10*i+5]+=cimag(dxiLM20dt);

			dydt[offset+10+10*i+6]+=creal(dxiLM21dt);
			dydt[offset+10+10*i+7]+=cimag(dxiLM21dt);

			dydt[offset+10+10*i+8]+=creal(dxiLM22dt);
			dydt[offset+10+10*i+9]+=cimag(dxiLM22dt);
		}
	}
#endif
	return GSL_SUCCESS;
}
