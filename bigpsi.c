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

int big_sc_time_evolution(double t,const double y[],double dydt[],void *data)
{
	struct bigpsi_t *bigpsi=data;
	struct params_t *params=bigpsi->params;
	struct configuration_t *config=params->config;

	double complex C;

	for(int c=0;c<(2+10*config->gridpoints)*bigpsi->nrpsis;c++)
		dydt[c]=0.0f;

	for(int c=0;c<bigpsi->nrpsis;c++)
	{
		int offset=c*(2+10*config->gridpoints);

		sc_time_evolution(t,&y[offset],&dydt[offset],&params[c]);
	}

	/*
		The level mixing due to the electric field
	*/

	if(config->laser==false)
		return GSL_SUCCESS;

	if((C=(2.0f/3.0f)*get_laser_intensity(config->milliwatts,config->duration,t,config))==0.0f)
		return GSL_SUCCESS;

	for(int d=0;d<bigpsi->nrpsis;d++)
	{
		int offset=d*(2+10*config->gridpoints);
		double complex gLM,dgLMdt;

		gLM=y[offset+0]+I*y[offset+1];

#define GETL(x) (params[x].L)
#define GETM() 	(params[0].M)

		dgLMdt=0.0f;
		for(int c=d-2;c<=d+2;c++)
		{
			int local_offset=c*(2+10*config->gridpoints);
			double complex gLprimeM;

			if((c<0)||(c>=bigpsi->nrpsis))
				continue;

			gLprimeM=y[local_offset+0]+I*y[local_offset+1];
			
			dgLMdt+=-I*C*Q(GETL(d),GETL(c),GETM(),0)*gLprimeM;
		}

		dgLMdt+=-I*(C/2.0f)*gLM;

		dydt[offset+0]+=creal(dgLMdt);
		dydt[offset+1]+=cimag(dgLMdt);

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

				if((c<=1)||(c>=bigpsi->nrpsis))
					continue;

				Lprime=GETL(c);

				phasediff=timephase((Lprime*(Lprime+1.0f)-L*(L+1.0f)),t,config);

				xiLprimeM2m2=y[local_offset+2+10*i]+I*y[local_offset+2+10*i+1];
				xiLprimeM2m1=y[local_offset+2+10*i+2]+I*y[local_offset+2+10*i+3];
				xiLprimeM20=y[local_offset+2+10*i+4]+I*y[local_offset+2+10*i+5];
				xiLprimeM21=y[local_offset+2+10*i+6]+I*y[local_offset+2+10*i+7];
				xiLprimeM22=y[local_offset+2+10*i+8]+I*y[local_offset+2+10*i+9];

#warning Controllare il segno dell'accoppiamento del laser

				dxiLM2m2dt+=-I*phasediff*C*Q(L,Lprime,GETM(),-2)*xiLprimeM2m2;
				dxiLM2m1dt+=-I*phasediff*C*Q(L,Lprime,GETM(),-1)*xiLprimeM2m1;
				dxiLM20dt+=-I*phasediff*C*Q(L,Lprime,GETM(),0)*xiLprimeM20;
				dxiLM21dt+=-I*phasediff*C*Q(L,Lprime,GETM(),1)*xiLprimeM21;
				dxiLM22dt+=-I*phasediff*C*Q(L,Lprime,GETM(),2)*xiLprimeM22;
			}

			xiLM2m2=y[offset+2+10*i]+I*y[offset+2+10*i+1];
			xiLM2m1=y[offset+2+10*i+2]+I*y[offset+2+10*i+3];
			xiLM20=y[offset+2+10*i+4]+I*y[offset+2+10*i+5];
			xiLM21=y[offset+2+10*i+6]+I*y[offset+2+10*i+7];
			xiLM22=y[offset+2+10*i+8]+I*y[offset+2+10*i+9];

#warning Controllare il segno dell'accoppiamento del laser

			dxiLM2m2dt+=-I*(C/2.0f)*xiLM2m2;
			dxiLM2m1dt+=-I*(C/2.0f)*xiLM2m1;
			dxiLM20dt+=-I*(C/2.0f)*xiLM20;
			dxiLM21dt+=-I*(C/2.0f)*xiLM21;
			dxiLM22dt+=-I*(C/2.0f)*xiLM22;

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
			psi->y[offset+2+10*c]=psi->y[offset+2+10*c+1]=0.0f;
			psi->y[offset+2+10*c+2]=psi->y[offset+2+10*c+3]=0.0f;
			psi->y[offset+2+10*c+4]=psi->y[offset+2+10*c+5]=0.0f;
			psi->y[offset+2+10*c+6]=psi->y[offset+2+10*c+7]=0.0f;
			psi->y[offset+2+10*c+8]=psi->y[offset+2+10*c+9]=0.0f;
		}
	}

	psi->t=config->starttime;

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

#warning Modify these parameters!

	psi->driver=gsl_odeiv2_driver_alloc_y_new(&psi->sys,gsl_odeiv2_step_rkf45,1e-9,1e-9,1e-9);

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

		gsl_odeiv2_driver_free(psi->driver);

		free(psi);
	}
}

#if 0
void bigpsi_serialize(struct bigpsi_t *psi,FILE *out)
{
	struct configuration_t *config=psi->config;

	double cutoff;
	int gridpoints;

	int c,ydim=(2+10*config->gridpoints)*psi->nrpsis;
	
	cutoff=config->cutoff;
	gridpoints=config->gridpoints;
	
	fwrite(&cutoff,sizeof(double),1,out);
	fwrite(&gridpoints,sizeof(int),1,out);

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

bool bigpsi_deserialize(FILE *in,struct bigpsi_t *psi,struct configuration_t *config)
{
	double cutoff;
	int gridpoints;

	int c,ydim=(2+10*config->gridpoints)*psi->nrpsis;
	
	fread(&cutoff,sizeof(double),1,in);
	fread(&gridpoints,sizeof(double),1,in);

	if((fabs(config->cutoff-cutoff)>=0.00001f)||(config->gridpoints!=gridpoints))
	{
		printf("Error: trying to load a saved state generated with gridpoints=%d, cutoff=%f\n",gridpoints,cutoff);
		printf("While the current code has: points=%d, gridcutoff=%f\n",config->gridpoints,config->cutoff);
		return false;
	}

	fread(&psi->nrpsis,sizeof(int),1,in);
	fread(&psi->dim,sizeof(int),1,in);

	bigpsi_init(psi,psi->nrpsis,0,0);

	fread(psi->y,sizeof(double),ydim,in);
	fread(&psi->t,sizeof(double),1,in);

	for(c=0;c<psi->nrpsis;c++)
	{
		fread(&psi->params[c].L,sizeof(int),1,in);
		fread(&psi->params[c].M,sizeof(int),1,in);
	}
	
	return true;
}
#endif

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

		num+=L*pow(norm_qp(psi->t,&psi->y[offset],config),2.0f);
		den+=pow(norm_qp(psi->t,&psi->y[offset],config),2.0f);
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

		ret+=pow(norm(psi->t,&psi->y[offset],config),2.0f);
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

		ret+=pow(norm_qp(psi->t,&psi->y[offset],config),2.0f);
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

		ret+=pow(norm_phonons(psi->t,&psi->y[offset],config),2.0f);
	}
	
	return sqrtf(ret);
}

void bigpsi_normalize(struct bigpsi_t *psi,double *normalization_error)
{
	struct configuration_t *config=psi->config;

	double norm=total_norm(psi);

	for(int c=0;c<(2+10*config->gridpoints)*psi->nrpsis;c++)
		psi->y[c]/=norm;
	
	if(normalization_error!=NULL)
		*normalization_error=norm;
}

void bigpsi_apply_step(struct bigpsi_t *psi,double ti,double *normalization_error)
{
	gsl_odeiv2_driver_apply(psi->driver,&psi->t,ti,psi->y);

	//bigpsi_normalize(psi,normalization_error);
}
