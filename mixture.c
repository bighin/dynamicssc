#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <sys/time.h>

#include "mixture.h"
#include "config.h"
#include "dsc.h"
#include "laser.h"
#include "bigpsi.h"
#include "cos.h"
#include "observables.h"

int do_run(int L,int M,struct info_t *info,bool silent,struct configuration_t *config)
{
	struct bigpsi_t *psi;
	int c,timedivs;

	char fname[1024];
	FILE *details;

	snprintf(fname,1024,"%s.singlestate.L%d.M%d.dat",config->prefix,L,M);
	fname[1023]='\0';

	if(!(details=fopen(fname,"w+")))
		fprintf(stderr,"Couldn't open %s for writing.\n",fname);

	fprintf(details,"# Running simulation with L=%d, M=%d, as part of a statistical mixture.\n",L,M);
	fprintf(details,"# <Time> <Intensity> <BathIntensity> <Re(L=0)> <Im(L=0)> ... <Re(L=Lmax)> <Im(L=Lmax)> <Re(cos3d)> <Im(cos3d)> <Re(cos2d)> <Im(cos2d)> <AOS> <TotalNorm> <TotalNormQP> <NormQP (L=0)> ... <NormQP (L=Lmax)> <NormTotal (L=0)> ... <NormTotal (L=Lmax)>\n");

	psi=bigpsi_init(config,BIGPSI_INIT_FROM_VALUES,L,M);

	timedivs=(config->endtime-config->starttime)/config->timestep;

	for(c=0;c<=timedivs;c++)
	{
		double ti=config->starttime+c*(config->endtime-config->starttime)/((double)(timedivs));
		double completion,elapsed_time;
		struct timeval starttime,endtime;

		gettimeofday(&starttime,NULL);

		bigpsi_apply_step(psi,ti,NULL,config);

		info[c].t=ti;

		info[c].totalnorm=total_norm(psi);
		info[c].totalnorm_qp=total_norm_qp(psi);
		info[c].intensity=get_laser_intensity(config->fluence,config->duration,ti,config);

                if(config->ramp)
                        info[c].bath_intensity=adiabatic_ramp(ti,config);
                else
                        info[c].bath_intensity=1.0f;

		info[c].cos2d=(config->cos2d==true)?(costheta2d(psi,config)):(0.0f);

		info[c].cossquared=costhetasquared(psi,config);
		info[c].aos=get_aos(psi);
		info[c].torque=fabs(torque(psi,L,M,config));

		info[c].j2=molecular_rotational_energy(psi,config);
		info[c].l2=total_rotational_energy(psi,config);
		info[c].lambda2=bosons_rotational_energy(psi,config);

		gettimeofday(&endtime,NULL);

		elapsed_time=(endtime.tv_sec-starttime.tv_sec)*1000.0;
		elapsed_time+=(endtime.tv_usec-starttime.tv_usec)/1000.0;
		elapsed_time/=1000;

		if(silent==false)
		{

#define REPSI_OFFSET(L)	(L*(10+10*config->gridpoints))
#define IMPSI_OFFSET(L)	(1+L*(10+10*config->gridpoints))

			double reL,imL;
		
			printf("%f %f %f ",ti,get_laser_intensity(config->fluence,config->duration,ti,config),info[c].bath_intensity);

			for(int L=0;L<=2;L++)
			{
				double complex z;

				reL=psi->y[REPSI_OFFSET(L)];
				imL=psi->y[IMPSI_OFFSET(L)];

				z=(reL+I*imL)*timephase(-L*(L+1.0f),psi->t,config);

				printf("%f %f ",creal(z),cimag(z));
			}

			printf("%f %f %f %f %f %f",creal(info[c].cos2d),cimag(info[c].cos2d),creal(info[c].cossquared),cimag(info[c].cossquared),info[c].totalnorm,info[c].totalnorm_qp);
		}
	
		completion=100.0f*(ti-config->starttime)/(config->endtime-config->starttime);

		printf("# L=%d, M=%d, norm=%f, localdensity=%f, %f%%, time=%f\n",L,M,info[c].totalnorm,config->density*info[c].bath_intensity,completion,elapsed_time);
		fprintf(details,"%f %f %f ",ti,get_laser_intensity(config->fluence,config->duration,ti,config),info[c].bath_intensity);

		{
			double reL,imL;
			double complex z;

			for(int L=0;L<config->maxl;L++)
			{
				reL=psi->y[REPSI_OFFSET(L)];
				imL=psi->y[IMPSI_OFFSET(L)];

				z=(reL+I*imL)*timephase(-L*(L+1.0f),psi->t,config);

				fprintf(details,"%f %f ",creal(z),cimag(z));
			}			
		}

		fprintf(details,"%f %f %f %f %f %f %f ",creal(info[c].cos2d),cimag(info[c].cos2d),creal(info[c].cossquared),cimag(info[c].cossquared),info[c].aos,info[c].totalnorm,info[c].totalnorm_qp);

		for(int d=0;d<=3;d++)
		{
			int offset=d*(10+10*config->gridpoints);
		
			fprintf(details,"%f ",norm_qp(ti,&(psi->y[offset]),&(psi->params[d]),config));
		}

		for(int d=0;d<config->maxl;d++)
		{
			int offset=d*(10+10*config->gridpoints);
		
			fprintf(details,"%f ",norm_qp(ti,&(psi->y[offset]),&(psi->params[d]),config));
		}

		for(int d=0;d<config->maxl;d++)
		{
			int offset=d*(10+10*config->gridpoints);
		
			fprintf(details,"%f ",norm(ti,&(psi->y[offset]),&(psi->params[d]),config));
		}

		fprintf(details,"\n");

		fflush(stdout);
		fflush(details);
	}

	bigpsi_fini(psi);

	if(details)
		fclose(details);

	return 0;
}

int compar_double(const void *a,const void *b)
{
	double first=*((double *)(a));
	double second=*((double *)(b));

	if(first>second)
		return 1;

	if(first<second)
		return -1;
	
	return 0;
}

double pctdiff(double a,double b)
{
	return 2.0f*fabs(a-b)/(a+b);
}

void debug_compare_weights(double *hardcoded_weights,int nr_states,struct configuration_t *config)
{
	double *calculated_weights=malloc(sizeof(double)*config->mixture_nr_states);

	for(int c=0;c<config->mixture_nr_states;c++)
		calculated_weights[c]=config->mixture_weights[c];

	qsort(hardcoded_weights,nr_states,sizeof(double),compar_double);
	qsort(calculated_weights,nr_states,sizeof(double),compar_double);

#define MAX(a,b)	(((a)>(b))?(a):(b))

	for(int c=0;c<MAX(nr_states,config->mixture_nr_states);c++)
	{
		printf("Hardcoded: ");

		if(c<nr_states)
			printf("%f, ",hardcoded_weights[c]);
		else
			printf("N/A ");

		printf("\tCalculated: ");

		if(c<config->mixture_nr_states)
			printf("%f ",calculated_weights[c]);
		else
			printf("N/A ");
	
		if((c<nr_states)&&(c<config->mixture_nr_states))
			printf("\tDiff: %f",pctdiff(hardcoded_weights[c],calculated_weights[c]));
		else
			printf("\tDiff: N/A");
		

		printf("\n");
	}

	exit(0);
}

void do_mixture(struct configuration_t *config)
{
	char outfile[1024];

	int timedivs;
	struct info_t **infos;
	int c,d;

	FILE *out;

	if(config->maxl<=6)
	{
		fprintf(stderr,"The 'maxl' parameter should be at least 6 for a statistical mixture. Exiting.\n");
		return;
	}

	if(config->mixture_nr_states==0)
	{
		fprintf(stderr,"The mixture weights have not been calculated for the molecule. Please check weather molecules.conf has been loaded and has the correct entry. Exiting.\n");
		return;
	}

	snprintf(outfile,1024,"%s.mixture.dat",config->prefix);
	outfile[1023]='\0';

	if(!(out=fopen(outfile,"w+")))
	{
		printf("Couldn't open %s for writing.\n",outfile);
		return;
	}

	infos=malloc(config->mixture_nr_states*sizeof(struct info_t *));

	timedivs=(config->endtime-config->starttime)/config->timestep;

	for(c=0;c<config->mixture_nr_states;c++)
		infos[c]=malloc(sizeof(struct info_t)*(1+timedivs));

#pragma omp parallel for

	for(c=0;c<config->mixture_nr_states;c++)
	{
		do_run(config->mixture_states[c][0],config->mixture_states[c][1],infos[c],false,config);
	}

	printf("# <Time> <Intensity (arbitrary units)> <BathIntensity> <Re(cos2d)> <Im(cos2d)> <Re(cos^2)> <Im(cos^2)> <AverageOccupiedState> <WeightedTotalNorm> <TotalOfWeights> <Torque> <J^2> <\\Lambda^2> <L^2>\n");
	fprintf(out,"# <Time> <Intensity (arbitrary units)> <BathIntensity> <Re(cos2d)> <Im(cos2d)> <Re(cos^2)> <Im(cos^2)> <AverageOccupiedState> <WeightedTotalNorm> <TotalOfWeights> <Torque> <J^2> <\\Lambda^2> <L^2>\n");

	for(d=0;d<=timedivs;d++)
	{
		double complex aceven=0.0f,cseven=0.0f;
		double complex acodd=0.0f,csodd=0.0f;
		double complex aosodd=0.0f,aoseven=0.0f;
		double trqodd=0.0,trqeven=0.0f;
		double j2odd=0.0,j2even=0.0f;
		double l2odd=0.0,l2even=0.0f;
		double lambda2odd=0.0,lambda2even=0.0f;

		double complex ac,cs;
		double aos;
		double trq;
		double j2,l2,lambda2;
		double atn=0.0f,tw=0.0f;

		double ea=config->mixture_even_abundance;
		double oa=config->mixture_odd_abundance;

		for(c=0;c<config->mixture_nr_states;c++)
		{
			int statel=config->mixture_states[c][0];

			if((statel%2)==0)
			{
				aceven+=config->mixture_weights[c]*infos[c][d].cos2d;
				cseven+=config->mixture_weights[c]*infos[c][d].cossquared;
				aoseven+=config->mixture_weights[c]*infos[c][d].aos;
				trqeven+=config->mixture_weights[c]*infos[c][d].torque;

				j2even+=config->mixture_weights[c]*infos[c][d].j2;
				l2even+=config->mixture_weights[c]*infos[c][d].l2;
				lambda2even+=config->mixture_weights[c]*infos[c][d].lambda2;
			}
			else
			{
				acodd+=config->mixture_weights[c]*infos[c][d].cos2d;
				csodd+=config->mixture_weights[c]*infos[c][d].cossquared;
				aosodd+=config->mixture_weights[c]*infos[c][d].aos;
				trqodd+=config->mixture_weights[c]*infos[c][d].torque;

				j2odd+=config->mixture_weights[c]*infos[c][d].j2;
				l2odd+=config->mixture_weights[c]*infos[c][d].l2;
				lambda2odd+=config->mixture_weights[c]*infos[c][d].lambda2;
			}

			atn+=infos[c][d].totalnorm;
			tw+=config->mixture_weights[c];
		}

		ac=(ea*aceven+oa*acodd)/(ea+oa);
		cs=(ea*cseven+oa*csodd)/(ea+oa);
		aos=(ea*aoseven+oa*aosodd)/(ea+oa);
		trq=(ea*trqeven+oa*trqodd)/(ea+oa);
		j2=(ea*j2even+oa*j2odd)/(ea+oa);
		l2=(ea*l2even+oa*l2odd)/(ea+oa);
		lambda2=(ea*lambda2even+oa*lambda2odd)/(ea+oa);

		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",infos[0][d].t,infos[0][d].intensity,infos[0][d].bath_intensity,creal(ac),cimag(ac),creal(cs),cimag(cs),aos,atn/config->mixture_nr_states,tw,trq,j2,lambda2,l2);
		fprintf(out,"%f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",infos[0][d].t,infos[0][d].intensity,infos[0][d].bath_intensity,creal(ac),cimag(ac),creal(cs),cimag(cs),aos,atn/config->mixture_nr_states,tw,trq,j2,lambda2,l2);
	}

	if(infos)
	{
		for(c=0;c<config->mixture_nr_states;c++)
			if(infos[c])
				free(infos[c]);
	
		free(infos);
	}

	if(out)
		fclose(out);

	printf("Done! The results have been written to %s.\n",outfile);
}
