#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <complex.h>

#include "dsc.h"
#include "laser.h"
#include "bigpsi.h"

double Q(int L,int Lprime,int M,int N)
{
	if((abs(N)>L)||(abs(N)>Lprime)||(abs(M)>L)||(abs(M)>Lprime))
		return 0.0f;
	
	return (sqrt(2.0f*L+1.0f)/sqrt(2.0f*Lprime+1.0f))*cg(2,0,L,M,Lprime,M)*cg(2,0,L,N,Lprime,N);
}

bool laser_is_on(struct configuration_t *config,double t)
{
	if(config->laser==false)
		return false;
	
	if(get_laser_intensity(config->milliwatts,config->duration,t,config)>1e-4)
		return true;
	
	return false;
}

int big_sc_time_evolution(double t,const double y[],double dydt[],void *data)
{
	struct bigpsi_t *bigpsi=data;
	struct params_t *params=bigpsi->params;
	struct configuration_t *config=params->config;

	double complex C;

	for(int c=0;c<(2+10*config->gridpoints)*bigpsi->nrpsis;c++)
		dydt[c]=0.0f;

	/*
		We apply the strong coupling evolution to each different L state...
	*/

	for(int c=0;c<bigpsi->nrpsis;c++)
	{
		int offset=c*(2+10*config->gridpoints);

		sc_time_evolution(t,&y[offset],&dydt[offset],&params[c]);
	}

	/*
		and then the level mixing due to the electric field!
	*/

	if(laser_is_on(config,t)==false)
		return GSL_SUCCESS;

	C=(2.0f/3.0f)*get_laser_intensity(config->milliwatts,config->duration,t,config);

#define GETL(x) (params[x].L)
#define GETM() 	(params[0].M)

	for(int d=0;d<bigpsi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);
		double complex gLM,dgLMdt;
		int L=GETL(d);

		gLM=(y[offset+0]+I*y[offset+1]);

		dgLMdt=0.0f;
		for(int c=d-2;c<=d+2;c++)
		{
			int local_offset=c*(2+10*config->gridpoints);
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

			dxiLM2m2dt=dxiLM2m1dt=dxiLM20dt=dxiLM21dt=dxiLM22dt=0.0f;
			for(int c=d-2;c<=d+2;c++)
			{
				int local_offset=c*(2+10*config->gridpoints);

				double complex xiLprimeM2m2,xiLprimeM2m1,xiLprimeM20,xiLprimeM21,xiLprimeM22;
				double complex phasediff;

				int Lprime;

				if((c<1)||(c>=bigpsi->nrpsis))
					continue;

				Lprime=GETL(c);

				phasediff=timephase((L*(L+1.0f)-Lprime*(Lprime+1.0f)),t,config);

				xiLprimeM2m2=y[local_offset+2+10*i]+I*y[local_offset+2+10*i+1];
				xiLprimeM2m1=y[local_offset+2+10*i+2]+I*y[local_offset+2+10*i+3];
				xiLprimeM20=y[local_offset+2+10*i+4]+I*y[local_offset+2+10*i+5];
				xiLprimeM21=y[local_offset+2+10*i+6]+I*y[local_offset+2+10*i+7];
				xiLprimeM22=y[local_offset+2+10*i+8]+I*y[local_offset+2+10*i+9];

				dxiLM2m2dt+=I*phasediff*C*Q(L,Lprime,GETM(),-2)*xiLprimeM2m2;
				dxiLM2m1dt+=I*phasediff*C*Q(L,Lprime,GETM(),-1)*xiLprimeM2m1;
				dxiLM20dt+=I*phasediff*C*Q(L,Lprime,GETM(),0)*xiLprimeM20;
				dxiLM21dt+=I*phasediff*C*Q(L,Lprime,GETM(),1)*xiLprimeM21;
				dxiLM22dt+=I*phasediff*C*Q(L,Lprime,GETM(),2)*xiLprimeM22;
			}

			xiLM2m2=y[offset+2+10*i]+I*y[offset+2+10*i+1];
			xiLM2m1=y[offset+2+10*i+2]+I*y[offset+2+10*i+3];
			xiLM20=y[offset+2+10*i+4]+I*y[offset+2+10*i+5];
			xiLM21=y[offset+2+10*i+6]+I*y[offset+2+10*i+7];
			xiLM22=y[offset+2+10*i+8]+I*y[offset+2+10*i+9];

			dxiLM2m2dt+=I*(C/2.0f)*xiLM2m2;
			dxiLM2m1dt+=I*(C/2.0f)*xiLM2m1;
			dxiLM20dt+=I*(C/2.0f)*xiLM20;
			dxiLM21dt+=I*(C/2.0f)*xiLM21;
			dxiLM22dt+=I*(C/2.0f)*xiLM22;

			dydt[offset+2+10*i]+=creal(dxiLM2m2dt);
			dydt[offset+2+10*i+1]+=cimag(dxiLM2m2dt);

			dydt[offset+2+10*i+2]+=creal(dxiLM2m1dt);
			dydt[offset+2+10*i+3]+=cimag(dxiLM2m1dt);

			dydt[offset+2+10*i+4]+=creal(dxiLM20dt);
			dydt[offset+2+10*i+5]+=cimag(dxiLM20dt);

			dydt[offset+2+10*i+6]+=creal(dxiLM21dt);
			dydt[offset+2+10*i+7]+=cimag(dxiLM21dt);

			dydt[offset+2+10*i+8]+=creal(dxiLM22dt);
			dydt[offset+2+10*i+9]+=cimag(dxiLM22dt);
		}
	}

	return GSL_SUCCESS;
}

struct bigpsi_t *bigpsi_init(struct configuration_t *config,int L,int M)
{
	struct bigpsi_t *psi;
	
	if(!(psi=malloc(sizeof(struct bigpsi_t))))
		return NULL;
	
	psi->nrpsis=config->maxl;
	psi->y=malloc((2+10*config->gridpoints)*psi->nrpsis*sizeof(double));
	psi->params=malloc(psi->nrpsis*sizeof(struct params_t));
	psi->config=config;
	psi->t=config->starttime;

	for(int d=0;d<psi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);

		if(d==L)
			psi->y[offset+0]=1.0f;
		else
			psi->y[offset+0]=0.0f;

		psi->y[offset+1]=0.0f;

		for(int c=0;c<config->gridpoints;c++)
		{
			psi->y[offset+2+10*c]=0.0f;
			psi->y[offset+2+10*c+1]=0.0f;
			psi->y[offset+2+10*c+2]=0.0f;
			psi->y[offset+2+10*c+3]=0.0f;
			psi->y[offset+2+10*c+4]=0.0f;
			psi->y[offset+2+10*c+5]=0.0f;
			psi->y[offset+2+10*c+6]=0.0f;
			psi->y[offset+2+10*c+7]=0.0f;
			psi->y[offset+2+10*c+8]=0.0f;
			psi->y[offset+2+10*c+9]=0.0f;
		}
	}

	for(int d=0;d<psi->nrpsis;d++)
	{
		psi->params[d].L=d;
		psi->params[d].M=M;
		psi->params[d].config=config;
	}

	psi->dim=psi->nrpsis*(2+10*config->gridpoints);

	psi->sys.function=big_sc_time_evolution;
	psi->sys.jacobian=NULL;
	psi->sys.dimension=psi->nrpsis*(2+10*config->gridpoints);
	psi->sys.params=psi;

	psi->driver=gsl_odeiv2_driver_alloc_y_new(&psi->sys,gsl_odeiv2_step_rkf45,config->hstart,config->epsabs,config->epsrel);

	psi->have_normalization_snapshot=false;
	psi->normalization_snapshot=malloc(sizeof(double)*psi->nrpsis);

	return psi;
}

void bigpsi_fini(struct bigpsi_t *psi)
{
	if(psi)
	{
		if(psi->y)
			free(psi->y);

		if(psi->params)
			free(psi->params);

		if(psi->normalization_snapshot)
			free(psi->normalization_snapshot);

		gsl_odeiv2_driver_free(psi->driver);

		free(psi);
	}
}

void bigpsi_serialize(struct bigpsi_t *psi,FILE *out)
{
	struct configuration_t *config=psi->config;

	int c,ydim=(2+10*config->gridpoints)*psi->nrpsis;

	fwrite(config,sizeof(struct configuration_t),1,out);

	fwrite(&psi->nrpsis,sizeof(int),1,out);
	fwrite(&psi->dim,sizeof(int),1,out);

	fwrite(psi->y,sizeof(double),ydim,out);
	fwrite(&psi->t,sizeof(double),1,out);

	for(c=0;c<psi->nrpsis;c++)
	{
		fwrite(&psi->params[c].L,sizeof(int),1,out);
		fwrite(&psi->params[c].M,sizeof(int),1,out);
	}
}

/*
	Here config must point to a (configuration_t *) struct, whose contents will be overwritten.
*/

struct bigpsi_t *bigpsi_deserialize(FILE *in,struct configuration_t *config)
{
	struct bigpsi_t *psi;
	int ydim;

	if(fread(config,sizeof(struct configuration_t),1,in)!=1)
		return NULL;

	ydim=(2+10*config->gridpoints)*config->maxl;
	psi=bigpsi_init(config,config->startl,config->startm);

	if(fread(&psi->nrpsis,sizeof(int),1,in)!=1)
		goto cleanup;
	
	if(fread(&psi->dim,sizeof(int),1,in)!=1)
		goto cleanup;

	if((psi->dim!=ydim)||(psi->nrpsis!=config->maxl))
	{
		fprintf(stderr,"Error loading a configuration from file! Data is inconsistent!\n");
		return NULL;
	}

	if(fread(psi->y,sizeof(double),ydim,in)!=ydim)
		return NULL;
	
	if(fread(&psi->t,sizeof(double),1,in)!=1)
		return NULL;

	for(int c=0;c<psi->nrpsis;c++)
	{
		int localL,localM;
		
		if(fread(&localL,sizeof(int),1,in)!=1)
			return NULL;
		
		if(fread(&localM,sizeof(int),1,in)!=1)
			return NULL;
	
		if((psi->params[c].L!=localL)||(psi->params[c].M!=localM))
		{
			fprintf(stderr,"Error loading a configuration from file! Data is inconsistent!\n");
			return NULL;
		}

		psi->params[c].config=config;
		psi->params[c].t=psi->t;
	}

	return psi;
	
	cleanup:
	
	bigpsi_fini(psi);
	return NULL;
}

double get_aos(struct bigpsi_t *psi)
{
	struct configuration_t *config=psi->config;

	int d;
	double num,den;
	
	num=den=0.0f;

	for(d=0;d<psi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);
		int L=psi->params[d].L;

		num+=L*pow(norm_qp(psi->t,&psi->y[offset],&psi->params[d],config),2.0f);
		den+=pow(norm_qp(psi->t,&psi->y[offset],&psi->params[d],config),2.0f);
	}
	
	return num/den;
}

double total_norm(struct bigpsi_t *psi)
{
	struct configuration_t *config=psi->config;

	int d;
	double ret=0.0f;

	for(d=0;d<psi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);

		ret+=pow(norm(psi->t,&psi->y[offset],&psi->params[d],config),2.0f);
	}
	
	return sqrtf(ret);
}

double total_norm_qp(struct bigpsi_t *psi)
{
	struct configuration_t *config=psi->config;

	int d;
	double ret=0.0f;

	for(d=0;d<psi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);

		ret+=pow(norm_qp(psi->t,&psi->y[offset],&psi->params[d],config),2.0f);
	}
	
	return sqrtf(ret);
}

double total_norm_phonons(struct bigpsi_t *psi)
{
	struct configuration_t *config=psi->config;

	int d;
	double ret=0.0f;

	for(d=0;d<psi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);

		ret+=pow(norm_phonons(psi->t,&psi->y[offset],&psi->params[d],config),2.0f);
	}
	
	return sqrtf(ret);
}

void bigpsi_normalize(struct bigpsi_t *psi,double *previousnorm)
{
	struct configuration_t *config=psi->config;

	/*
		It is better to normalize each L state in the wavefunction separately,
		due to better numerical accuracy.
	
		However, if the laser is on, or if we don't have a normalization snapshot,
		then we have 
	*/

	if((psi->have_normalization_snapshot==true)&&(laser_is_on(config,psi->t)==false))
	{
		double totalnorm=total_norm(psi);

		for(int c=0;c<psi->nrpsis;c++)
		{
			int offset=c*(2+10*config->gridpoints);
			double ratio;
			
			if(fabs(psi->normalization_snapshot[c])<=1e-10)
				continue;
			
			ratio=psi->normalization_snapshot[c]/norm(psi->t,&(psi->y[offset]),&psi->params[c],config);
			
			for(int d=0;d<(2+10*config->gridpoints);d++)
				psi->y[offset+d]*=ratio;
		}

		if(previousnorm!=NULL)
			*previousnorm=totalnorm;
	}
	else
	{
		double lnorm=total_norm(psi);

		for(int c=0;c<(2+10*config->gridpoints)*psi->nrpsis;c++)
			psi->y[c]/=lnorm;

		if(previousnorm!=NULL)
			*previousnorm=lnorm;

		/*
			We take a normalization snapshot, in case it could be used in the next stage.
		*/

		for(int c=0;c<psi->nrpsis;c++)
		{
			int offset=c*(2+10*config->gridpoints);

			psi->normalization_snapshot[c]=norm(psi->t,&(psi->y[offset]),&psi->params[c],config);
		}

		psi->have_normalization_snapshot=true;
	}
}

void bigpsi_apply_step(struct bigpsi_t *psi,double ti,double *previousnorm,struct configuration_t *config)
{
	gsl_odeiv2_driver_apply(psi->driver,&psi->t,ti,psi->y);

	if(config->normalize==true)
		bigpsi_normalize(psi,previousnorm);
}
