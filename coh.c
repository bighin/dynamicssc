#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>

#include "cubature/cubature.h"

#include "coh.h"
#include "dsc.h"
#include "bigpsi.h"
#include "observables.h"
#include "laser.h"

int fDcross_coherent(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct dint_container_t *container=(struct dint_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f;
	
	k=x[0];

	f=conj(get_point(container->intre1,k)+I*get_point(container->intim1,k));
	f*=get_point(container->intre2,k)+I*get_point(container->intim2,k);
	f/=fscale(k,container->L,config)*fscale(k,container->Lprime,config);

	switch(container->mode)
	{
		case DINT_MODE_PLAIN:
		break;

		case DINT_MODE_OMEGAK:
		f*=omegak(k,config);
		break;

		case DINT_MODE_VK:
		f*=pow(V2(k,container->localdensity,config)/W(k,config),2.0f);
		break;

		default:
		fprintf(stderr,"Unknown integration mode in Dcross()!\n");
		exit(0);
		break;
	}

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

/*
	Dcross_coherent() calculates an integral involving two phonon distribution functions as defined
	when using the coherent state Ansatz, more precisely:

	\sum_k \alpha^{* L}_{k 2 n} \alpha^{L'}_{k 2 n'} f(k)

	Note that the parameters L, Lprime, n and nprime are specified as arguments, Moreover one has

	f(k) = 1		if mode == DINT_MODE_PLAIN
	f(k) = omega(k)		if mode == DINT_MODE_OMEGAK
	f(k) = (V2(k)/W(k))^2	if mode == DINT_MODE_VK
*/

double complex Dcross_coherent(struct bigpsi_t *psi,int L,int Lprime,int n,int nprime,int mode,struct configuration_t *config)
{
	struct dint_container_t container;
	double *x,*y1re,*y1im,*y2re,*y2im;

	double xmin,xmax,res[2],err[2],localdensity;
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetL=L*(10+10*config->gridpoints);
	int offsetLprime=Lprime*(10+10*config->gridpoints);

	x=malloc(sizeof(double)*config->gridpoints);
	y1re=malloc(sizeof(double)*config->gridpoints);
	y1im=malloc(sizeof(double)*config->gridpoints);
	y2re=malloc(sizeof(double)*config->gridpoints);
	y2im=malloc(sizeof(double)*config->gridpoints);

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;
		int noffset,nprimeoffset;

		double complex phase1,phase2;

                x[c]=k;

		noffset=(n+2)*2;
		nprimeoffset=(nprime+2)*2;

		phase1=timephase(-(omegak(k,config)+6.0f),psi->t,config);
		phase2=timephase(-(omegak(k,config)+6.0f),psi->t,config);

                y1re[c]=creal(phase1*(psi->y[offsetL+10+10*c+noffset]+I*psi->y[offsetL+10+10*c+noffset+1]));
                y1im[c]=cimag(phase1*(psi->y[offsetL+10+10*c+noffset]+I*psi->y[offsetL+10+10*c+noffset+1]));

                y2re[c]=creal(phase2*(psi->y[offsetLprime+10+10*c+nprimeoffset]+I*psi->y[offsetLprime+10+10*c+nprimeoffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+10+10*c+nprimeoffset]+I*psi->y[offsetLprime+10+10*c+nprimeoffset+1]));
	}

        container.intre1=init_interpolation(x,y1re,config->gridpoints);
        container.intim1=init_interpolation(x,y1im,config->gridpoints);
        container.intre2=init_interpolation(x,y2re,config->gridpoints);
        container.intim2=init_interpolation(x,y2im,config->gridpoints);
	
	container.L=L;
	container.Lprime=Lprime;
	container.mode=mode;
	container.config=config;

        if(config->ramp==true)
                localdensity=config->density*adiabatic_ramp(psi->t,config);
        else
                localdensity=config->density;

	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fDcross_coherent,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intre1);
	fini_interpolation(container.intim1);
	fini_interpolation(container.intre2);
	fini_interpolation(container.intim2);

	if(x)	free(x);
	if(y1re) free(y1re);
	if(y1im) free(y1im);
	if(y2re) free(y2re);
	if(y2im) free(y2im);

	return res[0]+I*res[1];
}

int fDsingle_coherent(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        /* The context from which we read the global variables */

        struct dint_container_t *container=(struct dint_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables and other auxiliary variables */

        double k;
	double complex f;
	
	k=x[0];

	f=get_point(container->intre2,k)+I*get_point(container->intim2,k);
	f/=fscale(k,container->Lprime,config);

	switch(container->mode)
	{
		case DINT_MODE_PLAIN:
		case DINT_MODE_SUPERPLAIN:
		break;

		case DINT_MODE_OMEGAK:
		f*=omegak(k,config);
		break;

		case DINT_MODE_VK:
		f*=V2(k,container->localdensity,config)/W(k,config);
		break;

		case DINT_MODE_VK_OMEGAK:
		f*=omegak(k,config)*V2(k,container->localdensity,config)/W(k,config);
		break;

		case DINT_MODE_VK0:
		case DINT_MODE_SUPERPLAIN_VK0:
		f*=V2(k,container->localdensity,config);
		break;

		default:
		fprintf(stderr,"Unknown integration mode in Dcross()!\n");
		exit(0);
		break;
	}

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

/*
	Dsingle_coherent() calculates an integral involving one phonon distribution function as defined
	when using the coherent state Ansatz, more precisely:

	\sum_k \alpha^{L'}_{k 2 n} f(k)

	Note that the parameters L, Lprime and n are specified as arguments, Moreover one has

	f(k) = 1			if mode == DINT_MODE_SUPERPLAIN
	f(k) = V2(k)			if mode == DINT_MODE_SUPERPLAIN_VK0
*/

double complex Dsingle_coherent(struct bigpsi_t *psi,int L,int Lprime,int n,int mode,struct configuration_t *config)
{
	struct dint_container_t container;
	double *x,*y2re,*y2im;

	double xmin,xmax,res[2],err[2],localdensity;
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetLprime=Lprime*(10+10*config->gridpoints);

	x=malloc(sizeof(double)*config->gridpoints);
	y2re=malloc(sizeof(double)*config->gridpoints);
	y2im=malloc(sizeof(double)*config->gridpoints);

	if((mode!=DINT_MODE_SUPERPLAIN)&&(mode!=DINT_MODE_SUPERPLAIN_VK0))
	{
		fprintf(stderr,"Fatal error: wrong mode while calling Dsingle_coherent()\n");
		exit(0);
	}

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;
		int noffset;

		double complex phase2;

                x[c]=k;

		noffset=(n+2)*2;
		phase2=timephase(-(omegak(k,config)+6.0f),psi->t,config);

                y2re[c]=creal(phase2*(psi->y[offsetLprime+10+10*c+noffset]+I*psi->y[offsetLprime+10+10*c+noffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+10+10*c+noffset]+I*psi->y[offsetLprime+10+10*c+noffset+1]));
	}

        container.intre2=init_interpolation(x,y2re,config->gridpoints);
        container.intim2=init_interpolation(x,y2im,config->gridpoints);

	container.L=L;
	container.Lprime=Lprime;
	container.mode=mode;
	container.config=config;

        if(config->ramp==true)
                localdensity=config->density*adiabatic_ramp(psi->t,config);
        else
                localdensity=config->density;

	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fDsingle_coherent,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intre2);
	fini_interpolation(container.intim2);

	if(x)	free(x);
	if(y2re) free(y2re);
	if(y2im) free(y2im);

	return res[0]+I*res[1];
}

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

#warning Check the conventions of sigma

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(plus,mu,nu)*Dcross_coherent(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdaminus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int minus=1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(minus,mu,nu)*Dcross_coherent(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdazero(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int zero=1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(zero,mu,nu)*Dcross_coherent(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Vbeta(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	return creal(Dsingle_coherent(psi,L,L,0,DINT_MODE_SUPERPLAIN_VK0,config));
}

int sc_time_evolution_coherent(double t,const double y[],double dydt[],void *p)
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

	if((config->ramp==true)||(config->centrifugal))
	{
		printf("Some features (adiabatic ramp, finite temperature evolution, centrifugal distortion) are not supported at the moment when using coherent state evolution.\n");
		exit(0);
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

		dgndt+=-I*gn*Vbeta(psi,L,config)*timephase(L*(L+1.0f),t,config);
		dgndt+=I*gn*lambdasquared;

		dgndt+=-2.0f*I*gzero*lambdazero;
		dgndt+=-2.0f*I*gplus*lambdaminus;
		dgndt+=-2.0f*I*gminus*lambdaplus;

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
				dxi2mudt+=-I*V2(k,config->density,config)*timephase(omegak(k,config)+6.0f,t,config);

			dxi2mudt+=-2.0f*I*mu*gzero*xi2[2+mu];

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

int laser_time_evolution_coherent(double t,const double y[],double dydt[],void *data)
{
	struct bigpsi_t *bigpsi=data;
	struct params_t *params=bigpsi->params;
	struct configuration_t *config=params->config;

	double complex C;

	C=(2.0f/3.0f)*get_laser_intensity(config->fluence,config->duration,t,config);

#define GETL(x) (params[x].L)
#define GETM() 	(params[0].M)

#warning CHANGE FROM HERE!

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

	return GSL_SUCCESS;
}

double norm_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	double ret=0.0f;

	for(int n=-2;n<=2;n++)
	{
		double complex gn=y[2*(n+2)]+I*y[2*(n+2)+1];
	
		ret+=conj(gn)*gn;
	}

	return ret;
}

double norm_qp_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	return norm_coherent(t,y,params,config);
}

double norm_phonons_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
#warning PLEASE IMPLEMENT ME

	//Dcross_coherent(psi,Lprime,n,n,DINT_MODE_PLAIN,config);

	return 0.0f;
}

double complex fcosthetasquared_coherent(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f;
	double complex gL,gLprime,res;
	int n;

	int offsetL=L*(10+10*config->gridpoints);
	int offsetLprime=Lprime*(10+10*config->gridpoints);

	if(lambda==0)
		f=(1.0f/3.0f)*sqrt(4.0f*M_PI);

	if(lambda==2)
		f=(4.0/3.0)*sqrt(M_PI/5.0f);

	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);
	gLprime=(psi->y[offsetLprime+0]+I*psi->y[offsetLprime+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

	res=conj(gL)*gLprime*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
	res*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);

	for(n=-2;n<=2;n++)
	{
		double cgs;
		
		cgs=cg(L,M,lambda,0,Lprime,M)*cg(L,n,lambda,0,Lprime,n);
		
		if(fabs(cgs)>1e-10)
		{
			double complex tmp=1.0f;

			//tmp=Across(psi,L,Lprime,n,config)*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
			tmp*=cgs;
		
			res+=tmp;
		}
	}

	return f*res;
}

double complex costhetasquared_coh(struct bigpsi_t *psi,struct configuration_t *config)
{
	int L,Lprime,M,lambda;
	double complex total=0.0f;

	M=psi->params[0].M;

	for(L=0;L<config->maxl;L++)
	{
		for(Lprime=0;Lprime<config->maxl;Lprime++)
		{
			for(lambda=0;lambda<=2;lambda+=2)
			{
				total+=fcosthetasquared_coherent(psi,L,Lprime,M,lambda,config);
			}
		}
	}

	//if(config->altcos==true)
	//	return total/pow(total_norm_qp(psi),2.0f);

	return total;
}
