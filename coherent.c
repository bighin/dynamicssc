#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <gsl/gsl_errno.h>

#include "cubature/cubature.h"

#include "coherent.h"
#include "dsc.h"
#include "cos.h"
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

	double xmin,xmax,res[2],err[2];
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
                container.localdensity=config->density*adiabatic_ramp(psi->t,config);
        else
                container.localdensity=config->density;

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

#warning Check the conventions for sigma

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(plus,mu,nu)*Dcross_coherent(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdaminus(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int minus=-1;

	for(int mu=-2;mu<=2;mu++)
		for(int nu=-2;nu<=2;nu++)
			ret+=sigma_matrix(minus,mu,nu)*Dcross_coherent(psi,L,L,mu,nu,DINT_MODE_PLAIN,config);

	return ret;
}

double complex Lambdazero(struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	double complex ret=0.0f;
	int zero=0;

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

	double complex lambdaplus,lambdaminus,lambdazero,lambdasquared;
	double complex gplus,gminus,gzero;

	int L=params->L;

	gplus=Gplus(psi,L,config);
	gminus=Gminus(psi,L,config);
	gzero=Gzero(psi,L,config);

	lambdaplus=Lambdaplus(psi,L,config);
	lambdaminus=Lambdaminus(psi,L,config);
	lambdazero=Lambdazero(psi,L,config);

	lambdasquared=conj(lambdaplus)*lambdaplus+
		      conj(lambdaminus)*lambdaminus+
		      conj(lambdazero)*lambdazero;

	if((config->ramp==true)||(config->centrifugal==true))
	{
		printf("Some features (adiabatic ramp, centrifugal distortion) are not supported at the moment when using coherent state evolution.\n");
		exit(0);
	}

	for(int n=-2;n<=2;n++)
	{
		double complex gn,dgndt;
		int nplus,nminus;

		gn=y[2*(n+2)]+I*y[2*(n+2)+1];
		dgndt=0.0f;

		dgndt+=-I*gn*Vbeta(psi,L,config);
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

#warning The situation when the laser is on is much more complicated, one should derive the equations more carefully.

	for(int c=0;c<config->gridpoints;c++)
	{
		double gridstep=config->cutoff/config->gridpoints;
		double k=gridstep*c;

		double complex xi2[5];

		for(int mu=-2;mu<=2;mu++)
			xi2[2+mu]=y[10+10*c+(mu+2)*2]+I*y[10+10*c+(mu+2)*2+1];

		for(int mu=-2;mu<=2;mu++)
		{
			double complex dxi2mudt;
			int plus,minus,zero;

			dxi2mudt=0.0f;

			plus=+1;
			minus=-1;
			zero=0;

			if(mu==0)
				dxi2mudt+=-I*V2(k,config->density,config)*timephase(omegak(k,config)+6.0f,t,config);

			dxi2mudt+=2.0f*I*mu*gzero*xi2[2+mu];

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

			dydt[10+10*c+(mu+2)*2]=creal(dxi2mudt);
			dydt[10+10*c+(mu+2)*2+1]=cimag(dxi2mudt);	
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

	if(laser_is_on(config,t)==false)
		return GSL_SUCCESS;

	C=(2.0f/3.0f)*get_laser_intensity(config->fluence,config->duration,t,config);

#define GETL(x) (params[x].L)
#define GETM() 	(params[0].M)

	for(int d=0;d<bigpsi->nrpsis;d++)
	{
		int offset=d*(10+10*config->gridpoints);
		double complex gLMn,dgLMndt;
		int L=GETL(d);

		for(int n=-2;n<=2;n++)
		{
			gLMn=(y[offset+2*(n+2)]+I*y[offset+2*(n+2)+1]);

			dgLMndt=0.0f;
			for(int c=d-2;c<=d+2;c++)
			{
				int local_offset=c*(10+10*config->gridpoints);
				double complex gLprimeMn;
				int Lprime;

				if((c<0)||(c>=bigpsi->nrpsis))
					continue;

				Lprime=GETL(c);
				gLprimeMn=y[local_offset+2*(n+2)]+I*y[local_offset+2*(n+2)+1];

				dgLMndt+=I*C*Q(GETL(d),GETL(c),GETM(),0)*gLprimeMn*timephase((L*(L+1.0f)-Lprime*(Lprime+1.0f)),t,config);
			}

			dgLMndt+=I*(C/2.0f)*gLMn;

			dydt[offset+2*(n+2)]+=creal(dgLMndt);
			dydt[offset+2*(n+2)+1]+=cimag(dgLMndt);
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

	return sqrt(ret);
}

double norm_qp_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	return norm_coherent(t,y,params,config);
}

double norm_phonons_coherent(double t,const double y[],struct params_t *params,struct configuration_t *config)
{
	double complex total=0.0f;

	for(int n=-2;n<=2;n++)
	{
		struct dint_container_t container;
		double *x,*y1re,*y1im,*y2re,*y2im;

		double xmin,xmax,res[2],err[2],localdensity;
		int c;

		extern double relError;
		extern size_t maxEval;

		x=malloc(sizeof(double)*config->gridpoints);
		y1re=malloc(sizeof(double)*config->gridpoints);
		y1im=malloc(sizeof(double)*config->gridpoints);
		y2re=malloc(sizeof(double)*config->gridpoints);
		y2im=malloc(sizeof(double)*config->gridpoints);

		for(c=0;c<config->gridpoints;c++)
		{
			double gridstep=config->cutoff/config->gridpoints;
			double k=c*gridstep;
			int noffset;

			double complex phase;

			x[c]=k;

			noffset=(n+2)*2;
			phase=timephase(-(omegak(k,config)+6.0f),t,config);

			y1re[c]=creal(phase*(y[10+10*c+noffset]+I*y[10+10*c+noffset+1]));
			y1im[c]=cimag(phase*(y[10+10*c+noffset]+I*y[10+10*c+noffset+1]));

			y2re[c]=creal(phase*(y[10+10*c+noffset]+I*y[10+10*c+noffset+1]));
			y2im[c]=cimag(phase*(y[10+10*c+noffset]+I*y[10+10*c+noffset+1]));
		}

		container.intre1=init_interpolation(x,y1re,config->gridpoints);
		container.intim1=init_interpolation(x,y1im,config->gridpoints);
		container.intre2=init_interpolation(x,y2re,config->gridpoints);
		container.intim2=init_interpolation(x,y2im,config->gridpoints);

		container.L=container.Lprime=0;
		container.mode=DINT_MODE_PLAIN;
		container.config=config;

		if(config->ramp==true)
			localdensity=config->density*adiabatic_ramp(t,config);
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

		total+=res[0]+I*res[1];
	}

	return norm_coherent(t,y,params,config)*total;
}

double complex fcosthetasquared_coherent(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f;
	double complex res=0.0f;

	int offsetL=L*(10+10*config->gridpoints);
	int offsetLprime=Lprime*(10+10*config->gridpoints);

	if(lambda==0)
		f=1.0f/3.0f;

	if(lambda==2)
		f=2.0/3.0;

	for(int n=-2;n<=2;n++)
	{
		double complex gLn,gLprimen,localres;

		localres=0.0f;

		gLn=(psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1])*timephase(-L*(L+1),psi->t,config);
		gLprimen=(psi->y[offsetLprime+2*(n+2)]+I*psi->y[offsetLprime+2*(n+2)+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

		localres=conj(gLn)*gLprimen*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f));
		localres*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);

		res+=localres;
	}

	return f*res;
}

double complex costhetasquared_coherent(struct bigpsi_t *psi,struct configuration_t *config)
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

	return total;
}

double complex fcostheta2d_coherent(struct bigpsi_t *psi,int L,int Lprime,int M,int lambda,struct configuration_t *config)
{
	double f=flambda(lambda);
	double complex gL,gLprime,res;

	for(int n=-2;n<=2;n++)
	{
		int offsetL=L*(10+10*config->gridpoints);
		int offsetLprime=Lprime*(10+10*config->gridpoints);

		gL=(psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1])*timephase(-L*(L+1),psi->t,config);
		gLprime=(psi->y[offsetLprime+2*(n+2)]+I*psi->y[offsetLprime+2*(n+2)+1])*timephase(-Lprime*(Lprime+1),psi->t,config);

		res=conj(gL)*gLprime*sqrtf((2.0f*L+1.0f)/(2.0f*Lprime+1.0f))*sqrtf((2.0f*lambda+1.0f)/(4*M_PI));
		res*=cg(L,M,lambda,0,Lprime,M)*cg(L,0,lambda,0,Lprime,0);
	}

	return f*res;
}

double complex costheta2d_coherent(struct bigpsi_t *psi,struct configuration_t *config)
{
	int L,Lprime,M,lambda;
	double complex total=0.0f;

	M=psi->params[0].M;

	for(L=0;L<config->maxl;L++)
	{
		for(Lprime=0;Lprime<L;Lprime++)
		{
			for(lambda=0;lambda<config->maxl;lambda+=2)
			{
				total+=fcostheta2d_coherent(psi,L,Lprime,M,lambda,config)+
				       fcostheta2d_coherent(psi,Lprime,L,M,lambda,config);
			}
		}

		for(lambda=0;lambda<config->maxl;lambda++)
		{
			total+=fcostheta2d_coherent(psi,L,L,M,lambda,config);
		}
	}

	return total;
}

double phonon_overlap(double t,const double yt[],double s,const double ys[],struct params_t *params,struct configuration_t *config)
{
	double complex total=0.0f;

	for(int n=-2;n<=2;n++)
	{
		struct dint_container_t container;
		double *x,*y1re,*y1im,*y2re,*y2im;

		double xmin,xmax,res[2],err[2],localdensity;
		int c;

		extern double relError;
		extern size_t maxEval;

		x=malloc(sizeof(double)*config->gridpoints);
		y1re=malloc(sizeof(double)*config->gridpoints);
		y1im=malloc(sizeof(double)*config->gridpoints);
		y2re=malloc(sizeof(double)*config->gridpoints);
		y2im=malloc(sizeof(double)*config->gridpoints);

		for(c=0;c<config->gridpoints;c++)
		{
			double gridstep=config->cutoff/config->gridpoints;
			double k=c*gridstep;
			int noffset;

			double complex phases,phaset;

			x[c]=k;

			noffset=(n+2)*2;
			phases=timephase(-(omegak(k,config)+6.0f),s,config);
			phaset=timephase(-(omegak(k,config)+6.0f),t,config);

			y1re[c]=creal(phases*(ys[10+10*c+noffset]+I*ys[10+10*c+noffset+1]));
			y1im[c]=cimag(phases*(ys[10+10*c+noffset]+I*ys[10+10*c+noffset+1]));

			y2re[c]=creal(phaset*(yt[10+10*c+noffset]+I*yt[10+10*c+noffset+1]));
			y2im[c]=cimag(phaset*(yt[10+10*c+noffset]+I*yt[10+10*c+noffset+1]));
		}

		container.intre1=init_interpolation(x,y1re,config->gridpoints);
		container.intim1=init_interpolation(x,y1im,config->gridpoints);
		container.intre2=init_interpolation(x,y2re,config->gridpoints);
		container.intim2=init_interpolation(x,y2im,config->gridpoints);

		container.L=container.Lprime=0;
		container.mode=DINT_MODE_PLAIN;
		container.config=config;

		if(config->ramp==true)
			localdensity=config->density*adiabatic_ramp(t,config);
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

		total+=res[0]+I*res[1];
	}

	return total;
}

double complex overlapS_coherent(struct bigpsi_t *psi,double *y0,double t0,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
	{
		struct params_t *params=&psi->params[L];
		double complex norms,normt,overlapts,localres;
		int offsetL;

		offsetL=L*(10+10*config->gridpoints);
		localres=0.0f;

		for(int n=-2;n<=2;n++)
		{
			double complex gns,gnt;

			gns=(y0[offsetL+2*(n+2)]+I*y0[offsetL+2*(n+2)+1])*timephase(-L*(L+1.0f),t0,config);
			gnt=(psi->y[offsetL+2*(n+2)]+I*psi->y[offsetL+2*(n+2)+1])*timephase(-L*(L+1.0f),psi->t,config);
		
			localres+=conj(gns)*gnt;
		}

		norms=norm_phonons_coherent(t0,&y0[offsetL],params,config);
		normt=norm_phonons_coherent(psi->t,&psi->y[offsetL],params,config);
		overlapts=phonon_overlap(psi->t,&psi->y[offsetL],t0,&y0[offsetL],params,config);

		localres*=cexp(-0.5*norms-0.5*normt+overlapts);

		ret+=localres;
	}

	return ret;
}
