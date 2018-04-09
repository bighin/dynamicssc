#include <stdio.h>
#include <math.h>
#include <complex.h>

#include "cubature/cubature.h"

#include "auxx.h"
#include "bigpsi.h"
#include "config.h"
#include "observables.h"

int fDcross(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
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
	}

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

/*
	Dcross() calculates an integral involving two phonon distribution functions, more precisely:

	\sum_k \alpha^{* L'}_{k 2 n} \alpha^{L'}_{k 2 n'} f(k)

	Note that the parameters L, Lprime, n and nprime are specified as arguments, Moreover one has

	f(k) = 1		if mode == DINT_MODE_PLAIN
	f(k) = omega(k)		if mode == DINT_MODE_OMEGAK
	f(k) = (V2(k)/W(k))^2	if mode == DINT_MODE_VK
*/

double complex Dcross(struct bigpsi_t *psi,int L,int Lprime,int n,int nprime,int mode,struct configuration_t *config)
{
	struct dint_container_t container;
	double *x,*y1re,*y1im,*y2re,*y2im;

	double xmin,xmax,res[2],err[2],localdensity;
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

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

		switch(n)
		{
			case -2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=0;
			break;

			case -1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=2;
			break;

			case 0:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+6.0f),psi->t,config);
			noffset=4;
			break;

			case 1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=6;
			break;

			case 2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=8;
			break;
			
			default:
			phase1=0.0f;
			noffset=0;
			fprintf(stderr,"Bug in Dcross()\n");
			exit(0);
		}

		switch(nprime)
		{
			case -2:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			nprimeoffset=0;
			break;

			case -1:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			nprimeoffset=2;
			break;

			case 0:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+6.0f),psi->t,config);
			nprimeoffset=4;
			break;

			case 1:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			nprimeoffset=6;
			break;

			case 2:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			nprimeoffset=8;
			break;
			
			default:
			phase2=0.0f;
			nprimeoffset=0;
			fprintf(stderr,"Bug in Dcross()\n");
			exit(0);
		}

                y1re[c]=creal(phase1*(psi->y[offsetL+2+10*c+noffset]+I*psi->y[offsetL+2+10*c+noffset+1]));
                y1im[c]=cimag(phase1*(psi->y[offsetL+2+10*c+noffset]+I*psi->y[offsetL+2+10*c+noffset+1]));

                y2re[c]=creal(phase2*(psi->y[offsetLprime+2+10*c+nprimeoffset]+I*psi->y[offsetLprime+2+10*c+nprimeoffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+2+10*c+nprimeoffset]+I*psi->y[offsetLprime+2+10*c+nprimeoffset+1]));
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

	hcubature(2,fDcross,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

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

int fDsingle(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
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
		break;

		case DINT_MODE_OMEGAK:
		f*=omegak(k,config);
		break;

		case DINT_MODE_VK:
		f*=V2(k,container->localdensity,config)/W(k,config);
		break;
	}

	fval[0]=creal(f);
	fval[1]=cimag(f);

	return 0;
}

/*
	Dsingle() calculates an integral involving one phonon distribution function, more precisely:

	g^*_L \sum_k \alpha^{L'}_{k 2 n} f(k)

	Note that the parameters L, Lprime and n are specified as arguments, Moreover one has

	f(k) = 1		if mode == DINT_MODE_PLAIN
	f(k) = omega(k)		if mode == DINT_MODE_OMEGAK
	f(k) = V2(k)/W(k)	if mode == DINT_MODE_VK
*/

double complex Dsingle(struct bigpsi_t *psi,int L,int Lprime,int n,int mode,struct configuration_t *config)
{
	struct dint_container_t container;
	double *x,*y2re,*y2im;

	double xmin,xmax,res[2],err[2],localdensity;
	double complex g;
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

	x=malloc(sizeof(double)*config->gridpoints);
	y2re=malloc(sizeof(double)*config->gridpoints);
	y2im=malloc(sizeof(double)*config->gridpoints);

	g=timephase(-L*(L+1.0f),psi->t,config)*(psi->y[offsetL+0]+I*psi->y[offsetL+1]);

        for(c=0;c<config->gridpoints;c++)
	{
                double gridstep=config->cutoff/config->gridpoints;
                double k=c*gridstep;
		int noffset;

		double complex phase2;

                x[c]=k;

		switch(n)
		{
			case -2:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=0;
			break;

			case -1:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=2;
			break;

			case 0:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+6.0f),psi->t,config);
			noffset=4;
			break;

			case 1:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=6;
			break;

			case 2:
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=8;
			break;
			
			default:
			phase2=0.0f;
			noffset=0;
			fprintf(stderr,"Bug in Dsingle()\n");
			exit(0);
		}

                y2re[c]=creal(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
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

	hcubature(2,fDsingle,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

	fini_interpolation(container.intre2);
	fini_interpolation(container.intim2);

	if(x)	free(x);
	if(y2re) free(y2re);
	if(y2im) free(y2im);

	return conj(g)*(res[0]+I*res[1]);
}

/*
	Ecross() calculates an integral involving two phonon distribution functions
	at different times, more precisely:

	\sum_k \alpha^{* L'}_{k 2 n} (t=t0) \alpha^{L'}_{k 2 n} (t) f(k)

	Note that the parameters L, Lprime and n are specified as arguments.
*/

double complex Ecross(struct bigpsi_t *psi,int L,int Lprime,int n,double *y0,double t0,struct configuration_t *config)
{
	struct dint_container_t container;
	double *x,*y1re,*y1im,*y2re,*y2im;

	double xmin,xmax,res[2],err[2],localdensity;
	int c;

	extern double relError;
	extern size_t maxEval;

	int offsetL=L*(2+10*config->gridpoints);
	int offsetLprime=Lprime*(2+10*config->gridpoints);

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

		double complex phase1,phase2;

                x[c]=k;

		switch(n)
		{
			case -2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),t0,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=0;
			break;

			case -1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),t0,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=2;
			break;

			case 0:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+6.0f),t0,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+6.0f),psi->t,config);
			noffset=4;
			break;

			case 1:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),t0,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)+4.0f),psi->t,config);
			noffset=6;
			break;

			case 2:
			phase1=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),t0,config);
			phase2=timephase(-(Lprime*(Lprime+1.0f)+omegak(k,config)-2.0f),psi->t,config);
			noffset=8;
			break;
			
			default:
			phase1=phase2=0.0f;
			noffset=0;
			fprintf(stderr,"Bug in Ecross()\n");
			exit(0);
		}

                y1re[c]=creal(phase1*(y0[offsetL+2+10*c+noffset]+I*y0[offsetL+2+10*c+noffset+1]));
                y1im[c]=cimag(phase1*(y0[offsetL+2+10*c+noffset]+I*y0[offsetL+2+10*c+noffset+1]));

		y2re[c]=creal(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
                y2im[c]=cimag(phase2*(psi->y[offsetLprime+2+10*c+noffset]+I*psi->y[offsetLprime+2+10*c+noffset+1]));
	}

        container.intre1=init_interpolation(x,y1re,config->gridpoints);
        container.intim1=init_interpolation(x,y1im,config->gridpoints);
        container.intre2=init_interpolation(x,y2re,config->gridpoints);
        container.intim2=init_interpolation(x,y2im,config->gridpoints);
	
	container.L=L;
	container.Lprime=Lprime;
	container.mode=DINT_MODE_PLAIN;
	container.config=config;

        if(config->ramp==true)
                localdensity=config->density*adiabatic_ramp(psi->t,config);
        else
                localdensity=config->density;

	container.localdensity=localdensity;

	xmin=0.0f;
	xmax=config->cutoff;

	hcubature(2,fDcross,&container,1,&xmin,&xmax,maxEval,0,10*relError,ERROR_INDIVIDUAL,res,err);

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

/*
	See Misha's PRX, Eq. (D4).
*/

double eta_sigma(int L, int lambda,int n,int nu)
{
	if(n==nu)
		return n*n;

	if(n==(nu+1))
		return 0.5*sqrt(lambda*(lambda+1)-nu*(nu+1))*sqrt(L*(L+1)-nu*(nu+1));

	if(n==(nu-1))
		return 0.5*sqrt(lambda*(lambda+1)-nu*(nu-1))*sqrt(L*(L+1)-nu*(nu-1));

	return 0.0f;
}

int fD0(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        struct dint_container_t *container=(struct dint_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables */

        double k=x[0];

        fval[0]=pow(V2(k,container->localdensity,config)/W(k,config),2.0f);

        return 0;
}

double complex D0(double t,struct configuration_t *config)
{
	struct dint_container_t container;
        double xmin,xmax,res,err;

	extern double relError;
	extern size_t maxEval;

        xmin=0.0f;
        xmax=config->cutoff;

	container.config=config;

        if(config->ramp==true)
                container.localdensity=config->density*adiabatic_ramp(t,config);
        else
                container.localdensity=config->density;

        hcubature(1,fD0,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,&res,&err);

        return res;
}

double complex rotational_energy_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
#ifdef INCLUDE_EDEF
	double complex intD0=D0(psi->t,config);
#endif

	int offsetL=L*(2+10*config->gridpoints);

	double complex gL;
	double complex A1,A2,B2b,B3,C1,C2a,C2b,C3;
	
	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);

	A1=conj(gL)*gL*L*(L+1);
	A2=0.0f;

#ifdef INCLUDE_EDEF
	A2=conj(gL)*gL*6.0*intD0;
#endif

	B2b=B3=C1=C2a=C2b=C3=0.0f;

	if(L>=1)
	{
		B2b=-6.0*Dsingle(psi,L,L,0,DINT_MODE_VK,config);

		for(int n=-2;n<=2;n++)
			B3+=2.0f*eta_sigma(L,2,0,n)*Dsingle(psi,L,L,n,DINT_MODE_VK,config);

		for(int n=-2;n<=2;n++)
			C1+=L*(L+1)*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);

		for(int n=-2;n<=2;n++)
		{
			double complex a,b;

			if((n!=+1)&&(n!=-1))
				continue;

			a=Dsingle(psi,L,L,n,DINT_MODE_VK,config);
			b=0.5*6.0f*a*conj(a);

			C2a+=b+conj(b);
		}

		for(int n=-2;n<=2;n++)
		{
			C2b+=6.0f*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);
#ifdef INCLUDE_EDEF
			C2b+=6.0*intD0*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);
#endif
		}

		for(int n=-2;n<=2;n++)
			for(int nprime=-2;nprime<=2;nprime++)
				C3+=-2.0f*Dcross(psi,L,L,n,nprime,DINT_MODE_PLAIN,config)*eta_sigma(L,2,nprime,n);
	}

	return A1+A2+B2b+conj(B2b)+B3+conj(B3)+C1+C2a+C2b+C3;
}

double complex rotational_energy(struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
		ret+=rotational_energy_L(L,psi,config);
	
	return ret;
}

double complex rcr(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	return rotational_energy_L(L,psi,config)/(L*(L+1));
}

double complex overlapS(struct bigpsi_t *psi,double *y0,double t0,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
	{
		double complex gs,gt,phase1,phase2;
		int offsetL=L*(2+10*config->gridpoints);

		phase1=timephase(-L*(L+1.0f),t0,config);
		phase2=timephase(-L*(L+1.0f),psi->t,config);

		gs=phase1*(y0[offsetL+0]+I*y0[offsetL+1]);
		gt=phase2*(psi->y[offsetL+0]+I*psi->y[offsetL+1]);

		ret+=conj(gs)*gt;

		for(int n=-2;n<=2;n++)
			ret+=Ecross(psi,L,L,n,y0,t0,config);
	}

	return ret;
}
