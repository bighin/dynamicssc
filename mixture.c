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
	fprintf(details,"# <Time> <Intensity> <BathIntensity> <Re(L=0)> <Im(L=0)> ... <Re(L=Lmax)> <Im(L=Lmax)> <Re(cos3d)> <Im(cos3d)> <Re(cos2d)> <Im(cos2d)> <AOS> <TotalNorm> <TotalNormQP> <NormQP (L=0)> <NormQP (L=Lmax)>\n");

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

		gettimeofday(&endtime,NULL);

		elapsed_time=(endtime.tv_sec-starttime.tv_sec)*1000.0;
		elapsed_time+=(endtime.tv_usec-starttime.tv_usec)/1000.0;
		elapsed_time/=1000;

		if(silent==false)
		{

#define REPSI_OFFSET(L)	(L*(2+10*config->gridpoints))
#define IMPSI_OFFSET(L)	(1+L*(2+10*config->gridpoints))

			double reL,imL;
		
			printf("%f %f %f ",ti,get_laser_intensity(config->fluence,config->duration,ti,config),info[c].bath_intensity);

			for(int L=0;L<=2;L++)
			{
				reL=psi->y[REPSI_OFFSET(L)];
				imL=psi->y[IMPSI_OFFSET(L)];

				printf("%f %f ",reL,imL);
			}

			printf("%f %f %f %f %f %f",creal(info[c].cos2d),cimag(info[c].cos2d),creal(info[c].cossquared),cimag(info[c].cossquared),info[c].totalnorm,info[c].totalnorm_qp);
		}
	
		completion=100.0f*(ti-config->starttime)/(config->endtime-config->starttime);

		printf("# L=%d, M=%d, norm=%f, localdensity=%f, %f%%, time=%f\n",L,M,info[c].totalnorm,config->density*info[c].bath_intensity,completion,elapsed_time);
		fprintf(details,"%f %f %f ",ti,get_laser_intensity(config->fluence,config->duration,ti,config),info[c].bath_intensity);

		{
			double reL,imL;

			for(int L=0;L<=config->maxl;L++)
			{
				reL=psi->y[REPSI_OFFSET(L)];
				imL=psi->y[IMPSI_OFFSET(L)];

				fprintf(details,"%f %f ",reL,imL);
			}			
		}

		fprintf(details,"%f %f %f %f %f %f %f ",creal(info[c].cos2d),cimag(info[c].cos2d),creal(info[c].cossquared),cimag(info[c].cossquared),info[c].aos,info[c].totalnorm,info[c].totalnorm_qp);

		for(int d=0;d<=3;d++)
		{
			int offset=d*(2+10*config->gridpoints);
		
			fprintf(details,"%f ",norm_qp(ti,&(psi->y[offset]),&(psi->params[d]),config));
		}

		for(int d=0;d<config->maxl;d++)
		{
			int offset=d*(2+10*config->gridpoints);
		
			fprintf(details,"%f ",norm_qp(ti,&(psi->y[offset]),&(psi->params[d]),config));
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

void do_mixture(struct configuration_t *config)
{
	char outfile[1024];

	/*
		Statistical mixture for CS2, see EnsableAveragesCS2.nb and BW.nb for better precision.
	*/

#warning It would be nice to have an option to use the free weights (calculated without rotational constant renormalization)

#define USE_RENORMALIZED_WEIGHTS
#ifdef USE_RENORMALIZED_WEIGHTS

#define NR_CS2_STATES	(16)

	int cs2states[NR_CS2_STATES][2]=
	{
		{0, 0}, {2, 2}, {2, 1}, {2, 0}, {4, 4}, {4, 3}, {4, 2}, {4, 1},
		{4, 0}, {6, 6}, {6, 5}, {6, 4}, {6, 3}, {6, 2}, {6, 1}, {6, 0}
	};

  	double cs2weights[NR_CS2_STATES]=
	{
		0.237716, 0.226174, 0.226174, 0.113087, 0.0399576, 0.0399576,
		0.0399576, 0.0399576, 0.0199788, 0.00262157, 0.00262157, 0.00262157,
		0.00262157, 0.00262157, 0.00262157, 0.00131079
	};

	/*
		Statistical mixture for I2, see EnsableAveragesCS2.nb and BW.nb for better precision.
	*/

#define NR_I2_STATES	(28)

	int i2states[NR_I2_STATES][2]=
	{
		{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1},
		{3, 2}, {3, 3}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}, {5, 0},
		{5, 1}, {5, 2}, {5, 3}, {5, 4}, {5, 5}, {6, 0}, {6, 1}, {6, 2},
		{6, 3}, {6, 4}, {6, 5}, {6, 6}
	};

	double i2weights[NR_I2_STATES]=
	{
		0.165773, 0.142302, 0.284604, 0.0997231, 0.199446, 0.199446,
		0.0610023, 0.122005, 0.122005, 0.122005, 0.0304639, 0.0609277,
		0.0609277, 0.0609277, 0.0609277, 0.0132797, 0.0265595, 0.0265595,
		0.0265595, 0.0265595, 0.0265595, 0.00472587, 0.00945174, 0.00945174,
		0.00945174, 0.00945174, 0.00945174, 0.00945174
	};

	/*
		Statistical mixture for OCS.
	*/

#define NR_OCS_STATES	(15)

	int ocsstates[NR_OCS_STATES][2]=
	{
		{0, 0}, {1, 1}, {1, 0}, {2, 2}, {2, 1}, {2, 0}, {3, 3}, {3, 2},
		{3, 1}, {3, 0}, {4, 4}, {4, 3}, {4, 2}, {4, 1}, {4, 0}
	};

	double ocsweights[NR_OCS_STATES]=
	{
		0.25211, 0.290036, 0.145018, 0.0959652, 0.0959652, 0.0479826,
		0.0182645, 0.0182645, 0.0182645, 0.00913223, 0.00199955, 0.00199955,
		0.00199955, 0.00199955, 0.000999773
	};

#else

#define NR_CS2_STATES	(9)

	int cs2states[NR_CS2_STATES][2]=
	{
		{0, 0}, {2, 2}, {2, 1}, {2, 0}, {4, 4}, {4, 3}, {4, 2}, {4, 1}, {4, 0}
	};

  	double cs2weights[NR_CS2_STATES]=
	{
		0.702956, 0.11816, 0.11816, 0.0590798, 0.00036559,
		0.00036559,0.00036559, 0.00036559, 0.000182795
	};

#define NR_I2_STATES	(28)

	int i2states[NR_I2_STATES][2]=
	{
		{0, 0}, {1, 0}, {1, 1}, {2, 0}, {2, 1}, {2, 2}, {3, 0}, {3, 1},
		{3, 2}, {3, 3}, {4, 0}, {4, 1}, {4, 2}, {4, 3}, {4, 4}, {5, 0},
	        {5, 1}, {5, 2}, {5, 3}, {5, 4}, {5, 5}, {6, 0}, {6, 1}, {6, 2},
		{6, 3}, {6, 4}, {6, 5}, {6, 6}
	};

	double i2weights[NR_I2_STATES]=
	{
		0.269354, 0.203366, 0.406732, 0.115467, 0.230934, 0.230934,
		0.0495644, 0.0991289, 0.0991289, 0.0991289, 0.0159995, 0.031999,
		0.031999, 0.031999, 0.031999, 0.00390459, 0.00780918, 0.00780918,
		0.00780918, 0.00780918, 0.00780918, 0.000716585, 0.00143317,
		0.00143317, 0.00143317, 0.00143317, 0.00143317, 0.00143317
	};

#define NR_OCS_STATES	(10)

	int ocsstates[NR_OCS_STATES][2]=
	{
		{0, 0}, {1, 1}, {1, 0}, {2, 2}, {2, 1}, {2, 0},
		{3, 3}, {3, 2}, {3, 1}, {3, 0}
	};

	double ocsweights[NR_OCS_STATES]=
	{
		0.5895641146383894, 0.2537612526886281, 0.12688062634431405,
		0.011753138525620587, 0.011753138525620587, 0.005876569262810293,
		0.00011715118150224684, 0.00011715118150224684,0.00011715118150224684,
		0.00005857559075112342
	};

#endif

	double *weights;

	int nr_states,timedivs;
	struct info_t **infos;
	int c,d;

	FILE *out;

	if(config->maxl<6)
	{
		fprintf(stderr,"The 'maxl' parameter should be at least 6 for a statistical mixture. Exiting.\n");
		return;
	}

        switch(config->moleculetype)
        {
                case MOLECULE_I2:
		nr_states=NR_I2_STATES;
		weights=i2weights;
		break;

                case MOLECULE_CS2:
		nr_states=NR_CS2_STATES;
		weights=cs2weights;
                break;

                case MOLECULE_OCS:
		nr_states=NR_OCS_STATES;
		weights=ocsweights;
                break;

                default:
                fprintf(stderr,"Fatal error: unknown molecular species!\n");
                exit(0);
        }

	snprintf(outfile,1024,"%s.mixture.dat",config->prefix);
	outfile[1023]='\0';

	if(!(out=fopen(outfile,"w+")))
	{
		printf("Couldn't open %s for writing.\n",outfile);
		return;
	}

	infos=malloc(nr_states*sizeof(struct info_t *));

	timedivs=(config->endtime-config->starttime)/config->timestep;

	for(c=0;c<nr_states;c++)
		infos[c]=malloc(sizeof(struct info_t)*(1+timedivs));

#pragma omp parallel for

	for(c=0;c<nr_states;c++)
	{
	        switch(config->moleculetype)
	        {
	                case MOLECULE_I2:
			do_run(i2states[c][0],i2states[c][1],infos[c],false,config);

			break;

	                case MOLECULE_CS2:
			do_run(cs2states[c][0],cs2states[c][1],infos[c],false,config);
	                break;

	                case MOLECULE_OCS:
			do_run(ocsstates[c][0],ocsstates[c][1],infos[c],false,config);
	                break;
		}
	}

	printf("# <Time> <Intensity (arbitrary units)> <BathIntensity> <Re(cos2d)> <Im(cos2d)> <Re(cos^2)> <Im(cos^2)> <AverageOccupiedState> <WeightedTotalNorm> <TotalOfWeights> <Torque>\n");
	fprintf(out,"# <Time> <Intensity (arbitrary units)> <BathIntensity> <Re(cos2d)> <Im(cos2d)> <Re(cos^2)> <Im(cos^2)> <AverageOccupiedState> <WeightedTotalNorm> <TotalOfWeights> <Torque>\n");

	for(d=0;d<=timedivs;d++)
	{
		double complex aceven=0.0f,cseven=0.0f;
		double complex acodd=0.0f,csodd=0.0f;
		double complex aosodd=0.0f,aoseven=0.0f;
		double trqodd=0.0,trqeven=0.0f;
		double complex ac,cs;
		double aos;
		double trq;
		double atn=0.0f,tw=0.0f;

		for(c=0;c<nr_states;c++)
		{
			int statel;
			
		        switch(config->moleculetype)
		        {
		                case MOLECULE_I2:
				statel=i2states[c][0];
				break;

		                case MOLECULE_CS2:
				statel=cs2states[c][0];
		                break;

		                case MOLECULE_OCS:
				statel=ocsstates[c][0];
		                break;

				/*
					To silence silly gcc warnings.
				*/
	
				default:
				statel=0;
				break;
			}

			if((statel%2)==0)
			{
				aceven+=weights[c]*infos[c][d].cos2d;
				cseven+=weights[c]*infos[c][d].cossquared;
				aoseven+=weights[c]*infos[c][d].aos;
				trqeven+=weights[c]*infos[c][d].torque;
			}
			else
			{
				acodd+=weights[c]*infos[c][d].cos2d;
				csodd+=weights[c]*infos[c][d].cossquared;
				aosodd+=weights[c]*infos[c][d].aos;
				trqodd+=weights[c]*infos[c][d].torque;
			}

			atn+=infos[c][d].totalnorm;
			tw+=weights[c];
		}

	        switch(config->moleculetype)
	        {
	                case MOLECULE_I2:
			ac=(15.0f*aceven+21.0f*acodd)/(15.0f+21.0f);
			cs=(15.0f*cseven+21.0f*csodd)/(15.0f+21.0f);
			aos=(15.0f*aoseven+21.0f*aosodd)/(15.0f+21.0f);
			trq=(15.0f*trqeven+21.0f*trqodd)/(15.0f+21.0f);
			break;

	                case MOLECULE_CS2:
			ac=aceven;
			cs=cseven;
			aos=aoseven;
			trq=trqeven;
			break;

	                case MOLECULE_OCS:
			ac=(aceven+acodd)/(2.0f);
			cs=(cseven+csodd)/(2.0f);
			aos=(aoseven+aosodd)/(2.0f);
			trq=(trqeven+trqodd)/(2.0f);
			break;

			/*
				Again, just to silence silly gcc warnings.
			*/
		
			default:
			ac=cs=aos=trq=0.0f;
			break;
		}

		printf("%f %f %f %f %f %f %f %f %f %f %f\n",infos[0][d].t,infos[0][d].intensity,infos[0][d].bath_intensity,creal(ac),cimag(ac),creal(cs),cimag(cs),aos,atn/nr_states,tw,trq);
		fprintf(out,"%f %f %f %f %f %f %f %f %f %f %f\n",infos[0][d].t,infos[0][d].intensity,infos[0][d].bath_intensity,creal(ac),cimag(ac),creal(cs),cimag(cs),aos,atn/nr_states,tw,trq);
	}

	if(infos)
	{
		for(c=0;c<nr_states;c++)
			if(infos[c])
				free(infos[c]);
	
		free(infos);
	}

	if(out)
		fclose(out);

	printf("Done! The results have been written to %s.\n",outfile);
}
