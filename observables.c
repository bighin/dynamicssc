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
	Dcross() calculates an integral involving two phonon distribution functions, more precisely:

	\sum_k \alpha^{* L}_{k 2 n} \alpha^{L'}_{k 2 n'} f(k)

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
	Dsingle() calculates an integral involving one phonon distribution function, more precisely:

	g^*_L \sum_k \alpha^{L'}_{k 2 n} f(k)

	Note that the parameters L, Lprime and n are specified as arguments, Moreover one has

	f(k) = 1			if mode == DINT_MODE_PLAIN
	f(k) = omega(k)			if mode == DINT_MODE_OMEGAK
	f(k) = V2(k)/W(k)		if mode == DINT_MODE_VK
	f(k) = V2(k)			if mode == DINT_MODE_VK0
	f(k) = omega(k)*V2(k)/W(k)	if mode == DINT_MODE_VK_OMEGAK

	If mode == DINT_MODE_SUPERPLAIN then the integral calculated is

	\sum_k \alpha^{L'}_{k 2 n}
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

	if(mode==DINT_MODE_SUPERPLAIN)
		return res[0]+I*res[1];

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

int fD1(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        struct dint_container_t *container=(struct dint_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables */

        double k=x[0];

        fval[0]=omegak(k,config)*pow(V2(k,container->localdensity,config)/W(k,config),2.0f);

        return 0;
}

double complex D1(double t,struct configuration_t *config)
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

        hcubature(1,fD1,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,&res,&err);

        return res;
}

int fD2(unsigned ndim,const double *x,void *fdata,unsigned fdim,double *fval)
{
        struct dint_container_t *container=(struct dint_container_t *)(fdata);
	struct configuration_t *config=container->config;

        /* The integration variables */

        double k=x[0];

        fval[0]=pow(V2(k,container->localdensity,config),2.0f)/W(k,config);

        return 0;
}

double complex D2(double t,struct configuration_t *config)
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

        hcubature(1,fD2,&container,1,&xmin,&xmax,maxEval,0,relError,ERROR_INDIVIDUAL,&res,&err);

        return res;
}

double complex bosons_rotational_energy_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	int offsetL=L*(2+10*config->gridpoints);

	double complex gL;
	double norm2L;

	/*
		A2C2p1 contains A2 and the first part of C2
		C2p2 contains the second part of C2

		...and so on!
	*/

	double complex A2C2p1,C2p2,B2,C2p3;
	
	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);

	A2C2p1=C2p2=B2=C2p3=0.0f;

	norm2L=pow(norm(psi->t,&psi->y[offsetL],&psi->params[L],config),2.0f);

	A2C2p1=6.0*D0(psi->t,config)*norm2L;
	C2p2=6.0*(norm2L-conj(gL)*gL);

	if(L>=1)
	{
		for(int n=-2;n<=2;n++)
		{
			double complex B=0.0;
			
			if((n!=+1)&&(n!=-1))
				continue;

			if(fabs(gL)>1e-5)
				B=Dsingle(psi,L,L,n,DINT_MODE_VK,config)/conj(gL);

			C2p3+=6.0*B*conj(B);
		}

		B2+=-6.0*Dsingle(psi,L,L,0,DINT_MODE_VK,config);
	}

	return A2C2p1+C2p2+B2+conj(B2)+C2p3;
}

double complex bosons_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
		ret+=bosons_rotational_energy_L(L,psi,config);

	return ret;
}

double complex molecular_rotational_energy_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	int offsetL=L*(2+10*config->gridpoints);

	double norm2L;
	double complex A1C1,B3,C3;

	norm2L=pow(norm(psi->t,&psi->y[offsetL],&psi->params[L],config),2.0f);

	A1C1=B3=C3=0.0f;

	A1C1=norm2L*L*(L+1);

	if(L>=1)
	{
		for(int n=-2;n<=2;n++)
			for(int nprime=-2;nprime<=2;nprime++)
				C3+=-2.0f*eta_sigma(L,2,nprime,n)*Dcross(psi,L,L,nprime,n,DINT_MODE_PLAIN,config);

		for(int n=-2;n<=2;n++)
			B3+=2.0f*eta_sigma(L,2,0,n)*Dsingle(psi,L,L,n,DINT_MODE_VK,config);
	}

	return bosons_rotational_energy_L(L,psi,config)+A1C1+B3+conj(B3)+C3;
}

/*
	<J^2>_{lab} = <(J'-\Lambda)^2>_{rot}
*/

double complex molecular_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
		ret+=molecular_rotational_energy_L(L,psi,config);
	
	return ret;
}

/*
	<L^2>_{lab} = <L^2>_{rot}
*/

double complex total_rotational_energy(struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
	{
		int offsetL;
		double norm2L;

		offsetL=L*(2+10*config->gridpoints);
		norm2L=pow(norm(psi->t,&psi->y[offsetL],&psi->params[L],config),2.0f);

		ret+=norm2L*L*(L+1);
	}
	
	return ret;
}

double complex total_energy_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	int offsetL;
	double complex gL;
	double complex A4,A5,B4,B5,C4,C5;
	
	offsetL=L*(2+10*config->gridpoints);
	gL=(psi->y[offsetL+0]+I*psi->y[offsetL+1])*timephase(-L*(L+1),psi->t,config);

	A4=A5=B4=B5=C4=C5=0.0f;

	A4=D1(psi->t,config)*conj(gL)*gL;
	A5=-2.0f*D2(psi->t,config)*conj(gL)*gL;

	if(L>=1)
	{
		B4+=-Dsingle(psi,L,L,0,DINT_MODE_VK_OMEGAK,config);
		B5+=Dsingle(psi,L,L,0,DINT_MODE_VK0,config);

		for(int n=-2;n<=2;n++)
		{
			C4+=Dcross(psi,L,L,n,n,DINT_MODE_OMEGAK,config);
			C4+=D1(psi->t,config)*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);

			C5+=-2.0f*D2(psi->t,config)*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);
		}
	}

	return A4+A5+B4+conj(B4)+B5+conj(B5)+C4+C5;
}

double complex total_energy(struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex ret=0.0f;

	for(int L=0;L<config->maxl;L++)
	{
		ret+=molecular_rotational_energy_L(L,psi,config);
		ret+=total_energy_L(L,psi,config);
	}
	
	return ret;
}

double sigma_matrix(int i,int n,int nprime)
{
	if((abs(n)>2)||(abs(nprime)>2))
	{
		printf("Error in sigma_matrix(): wrong arguments.");
		exit(0);
	}

	switch(i)
	{
		case 1:

		if((n==-2)&&(nprime==-1))
			return -sqrtf(2.0f);

		if((n==-1)&&(nprime==0))
			return -sqrtf(3.0f);

		if((n==0)&&(nprime==1))
			return -sqrtf(3.0f);

		if((n==1)&&(nprime==2))
			return -sqrtf(2.0f);

		break;

		case 0:

		if(n==nprime)
			return -n;

		break;

		case -1:

		if((n==-1)&&(nprime==-2))
			return sqrtf(2.0f);
		
		if((n==0)&&(nprime==-1))
			return sqrtf(3.0f);

		if((n==1)&&(nprime==0))
			return sqrtf(3.0f);

		if((n==2)&&(nprime==1))
			return sqrtf(2.0f);

		break;
	}
	
	return 0.0f;
}

void debug_sigma_matrices(void)
{
	for(int i=-1;i<=1;i++)
	{
		for(int n=-2;n<=2;n++)
		{
			for(int nprime=-2;nprime<=2;nprime++)
			{
				printf("%f ",sigma_matrix(i,n,nprime));
			}

			printf("\n");
		}
		
		printf("\n\n");
	}
}

double complex JdotLambda_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex X3,Y3;
	int M;
	bool debugoutput,printout=false;
	
	X3=Y3=0.0f;
	M=psi->params[L].M;

	if(L>=1)
	{
		/*
			Enable the following option to have a verbose, debug output.
		*/
		
		debugoutput=false;

		for(int n=-2;n<=2;n++)
			if(fabs(Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config))>1e-8)
				printout=true;
		
		/*
			Calculation of X3: let's iterate over all n's
		*/

		for(int i=-1;i<=1;i++)
		{
			if(i!=0)
				continue;
			
			for(int n=-2;n<=2;n++)
			{
				double complex localres;
			
				/*
					The next line calculates the integral
			
					\sum_k \alpha^{* L}_{k 2 n} \alpha^{L}_{k 2 n}
			
					and saves it into the variable 'localres'.
				*/
			
				localres=Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);
			
				/*
					Here we multiply (notice the *= operator) the variable 'localres'
					by a factor
				
						|C_{L(n+i) 1-i}^{Ln}|^2
				*/
				
				localres*=pow(cg(L,n+i,1,-i,L,n),2.0f);
			
				/*
					We multiply the variable 'localres' by n
				*/
			
				localres*=n;
		
				/*
					Finally we multiply the variable 'localres' by sqrt(L(L+1)) * C_{LM 10}^{LM}
				*/
		
				localres*=sqrtf(L*(L+1.0f))*cg(L,M,1,0,L,M);
		
				/*
					We save everything in X3...
				*/
		
				X3+=localres;
		
				/*
					...and we print the result if needed.
				*/
		
				if((debugoutput==true)&&(printout==true))
					printf("X3(L=%d,i=%d,n=%d) = %f\n",L,i,n,creal(localres));
			}
		}

		/*
			Calculation of Y3: let's iterate over all i's and n's
		*/

		for(int i=-1;i<=1;i++)
		{
			for(int n=-2;n<=2;n++)
			{
				double complex localres;
				
				if(abs(n+i)>2)
					continue;
			
				/*
					The next line calculates the integral
			
					\sum_k \alpha^{* L}_{k 2 n} \alpha^{L}_{k 2 (n-i)}
			
					and saves it into the variable 'localres'.
				*/
	
				localres=Dcross(psi,L,L,n,n+i,DINT_MODE_PLAIN,config);
					
				/*
					Then we multiply by (-1)^i * C_{Ln 1-i}^{L (n-1)}
				*/

				localres*=cg(L,n,1,i,L,n+i);

				/*
					Then we multiply by (\sigma_{n (n-i)})_i * (n-i)
				*/

				localres*=sigma_matrix(i,n+i,n)*(-1)*(n+i);

				/*
					Finally we multiply the variable 'localres' by sqrt(L(L+1)) * C_{LM 10}^{LM}
				*/

				localres*=cg(L,M,1,0,L,M);

				/*
					We save everything in Y3...
				*/
		
				Y3+=localres;
		
				/*
					...and we print the result if needed.
				*/				

				if((debugoutput==true)&&(printout==true))
					printf("Y3(L=%d,i=%d,n=%d) = %f\n",L,i,n,creal(localres));
			}
		}
	}

	return X3-Y3;
}

double complex JdotLambda(struct bigpsi_t *psi,struct configuration_t *config)
{
        double complex ret=0.0f;

        for(int L=0;L<config->maxl;L++)
                ret+=JdotLambda_L(L,psi,config);

        return ret;
}

double complex Lambdaz2_rot_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	//int offsetL;
	double complex ret=0.0f;
	//double norm2L;

	//offsetL=L*(2+10*config->gridpoints);
	//norm2L=pow(norm(psi->t,&psi->y[offsetL],&psi->params[L],config),2.0f);

	for(int n=-2;n<=2;n++)
		ret+=n*n*Dcross(psi,L,L,n,n,DINT_MODE_PLAIN,config);
	
	return ret;
}

/*
	\langle \Lambda_z^2 \rangle_{rot}
*/

double complex Lambdaz2_rot(struct bigpsi_t *psi,struct configuration_t *config)
{
        double complex ret=0.0f;

        for(int L=0;L<config->maxl;L++)
                ret+=Lambdaz2_rot_L(L,psi,config);

        return ret;
}

double complex fA_L(int L,int n,int nprimeprime,int mu,struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex X3,Y3;
	int M;
	
	X3=Y3=0.0f;
	M=psi->params[L].M;

	if(L>=1)
	{
		for(int i=-1;i<=1;i++)
		{
			double complex localres=0.0f;

			//if(abs(n+i)>2)
			//	continue;

			if(n==mu)
			{
				localres=sqrtf(L*(L+1.0f))*pow(-1.0f,i+1);
				localres*=cg(L,n+i,1,-i,L,n);
				localres*=cg(L,M,1,0,L,M);
				localres*=n;
				localres*=cg(L,n+i,1,i,L,nprimeprime);
			}

			X3+=localres;
		}

		for(int i=-1;i<=1;i++)
		{
			double complex localres;
				
			//if(abs(n+i)>2)
			//	continue;

			localres=cg(L,M,1,0,L,M);
			localres*=cg(L,n,1,i,L,nprimeprime);
			localres*=-1.0f;
			localres*=sigma_matrix(i,mu,n);
			localres*=mu;

			Y3+=localres;
		}
	}

	return X3-Y3;
}

double complex Delta_JdotLambda_L(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex A[5][5][5];
	double complex intalpha[5][5];
	double complex ret=0.0f;

	for(int n=-2;n<=2;n++)
		for(int nprimeprime=-2;nprimeprime<=2;nprimeprime++)
			for(int mu=-2;mu<=2;mu++)
				A[2+n][2+nprimeprime][2+mu]=fA_L(L,n,nprimeprime,mu,psi,config);

	for(int n=-2;n<=2;n++)
		for(int nprime=-2;nprime<=2;nprime++)
			intalpha[2+n][2+nprime]=Dcross(psi,L,L,n,nprime,DINT_MODE_PLAIN,config);

	for(int n=-2;n<=2;n++)
		for(int nprime=-2;nprime<=2;nprime++)
			for(int nprimeprime=-2;nprimeprime<=2;nprimeprime++)
				for(int mu=-2;mu<=2;mu++)
					ret+=intalpha[2+n][2+nprime]*A[2+n][2+nprimeprime][2+mu]*conj(A[2+nprime][2+nprimeprime][2+mu]);

	return ret;
}

double complex Delta_JdotLambda(struct bigpsi_t *psi,struct configuration_t *config)
{
        double complex ret=0.0f;

        for(int L=0;L<config->maxl;L++)
                ret+=Delta_JdotLambda_L(L,psi,config);

        return ret;
}

double complex Jz_lab_L(int L,int M,struct bigpsi_t *psi,struct configuration_t *config)
{
	double complex A,B;
	
	A=B=0.0f;

	for(int i=-1;i<=1;i++)
	{
		for(int n=-2;n<=2;n++)
		{
			double cf=cg(L,M,1,0,L,M)*cg(L,0,1,-i,L,n)*sigma_matrix(i,n,0);
			
			if(fabs(cf)<1e-8)
				continue;

			A+=cf*Dsingle(psi,L,L,n,DINT_MODE_VK,config);
		}
	}

	for(int i=-1;i<=1;i++)
	{
		for(int n=-2;n<=2;n++)
		{
			for(int nprime=-2;nprime<=2;nprime++)
			{
				double cf=cg(L,M,1,0,L,M)*cg(L,nprime,1,-i,L,n)*sigma_matrix(i,n,nprime);
			
				if(fabs(cf)<1e-8)
					continue;
				
				B-=cf*Dcross(psi,L,L,n,nprime,DINT_MODE_PLAIN,config);
			}
		}
	}


	return A+conj(A)+B;
}

double complex Jz_lab(struct bigpsi_t *psi,struct configuration_t *config)
{
        double complex ret=0.0f;
	int M=psi->params[0].M;

        for(int L=0;L<config->maxl;L++)
                ret+=Jz_lab_L(L,M,psi,config);

        return M+ret;
}

double complex rcr(int L,struct bigpsi_t *psi,struct configuration_t *config)
{
	return molecular_rotational_energy_L(L,psi,config)/(L*(L+1));
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

#warning Only the molecular part of the overlap is here

		//for(int n=-2;n<=2;n++)
		//	ret+=Ecross(psi,L,L,n,y0,t0,config);
	}

	return ret;
}

double torque(struct bigpsi_t *psi,int L,int M,struct configuration_t *config)
{
	//double complex gLM;
	double cg;
	double complex integral;
	//int offsetL;

	if(L==0)
		return 0.0f;

	//offsetL=L*(2+10*config->gridpoints);
	//gLM=timephase(-L*(L+1.0f),psi->t,config)*(psi->y[offsetL+0]+I*psi->y[offsetL+1]);

	cg=1.0f*M/sqrtf(L*L+L);
	integral=conj(Dsingle(psi,L,L,1,DINT_MODE_VK0,config)+Dsingle(psi,L,L,-1,DINT_MODE_VK0,config));

	/*
		Note: we could also return the torque normalized by |g_{LM}|^2
	*/

	//return 2.0f*sqrtf(6.0f)*cg*cimag(integral)/(conj(gLM)*gLM);
	return 2.0f*sqrtf(6.0f)*cg*cimag(integral);
}
