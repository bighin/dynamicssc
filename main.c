#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/time.h>

#include "inih/ini.h"

#include "auxx.h"
#include "config.h"
#include "laser.h"
#include "bigpsi.h"
#include "cos.h"
#include "mixture.h"
#include "observables.h"
#include "molecules.h"

void print_header_g(FILE *out,int L,struct configuration_t *config)
{
	fprintf(out,"# Strong coupling time evolution\n");
	fprintf(out,"#\n");
	fprintf(out,"# L=%d\n",L);
	fprintf(out,"#\n");
	fprintf(out,"# Other parameters are:\n");
	fprintf(out,"# Cutoff: %f\n",(double)(config->cutoff));
	fprintf(out,"# Gridpoints: %i\n",config->gridpoints);

	if(!config->shapefile)
	{
		fprintf(out,"# Fluence (J/cm^2): %f\n",config->fluence);
		fprintf(out,"# Duration (ps): %f\n",config->duration);
	}
	else
	{
		fprintf(out,"# Pulse shape loaded from file '%s'\n",config->shapefile);
	}

	fprintf(out,"# Minimum time: %f\n",config->starttime);
	fprintf(out,"# Maximum time: %f\n",config->endtime);
	fprintf(out,"# Maximum L: %i\n",config->maxl);
	fprintf(out,"# M: %i\n",config->startm);
	fprintf(out,"#\n");
	fprintf(out,"# t Re(g) Im(g) Norm\n");
}

void print_header_alpha(FILE *out,int L,struct configuration_t *config)
{
	fprintf(out,"# Strong coupling time evolution\n");
	fprintf(out,"#\n");
	fprintf(out,"# L=%d\n",L);
	fprintf(out,"#\n");
	fprintf(out,"# Other parameters are:\n");
	fprintf(out,"# Cutoff: %f\n",((double)(config->cutoff)));
	fprintf(out,"# Gridpoints: %i\n",config->gridpoints);

	if(!config->shapefile)
	{
		fprintf(out,"# Fluence (J/cm^2): %f\n",config->fluence);
		fprintf(out,"# Duration (ps): %f\n",config->duration);
	}
	else
	{
		fprintf(out,"# Pulse shape loaded from file '%s'\n",config->shapefile);
	}

	fprintf(out,"# Minimum time: %f\n",config->starttime);
	fprintf(out,"# Maximum time: %f\n",config->endtime);
	fprintf(out,"# Maximum L: %i\n",config->maxl);
	fprintf(out,"# M: %i\n",config->startm);
	fprintf(out,"#\n");
	fprintf(out,"# t k Re(alpha2,-2) Im(alpha2,-2) Re(alpha2,-1) Im(alpha2,-1) Re(alpha2,0) Im(alpha2,0) Re(alpha2,1) Im(alpha2,1) Re(alpha2,2) Im(alpha2,2)\n");
}

void print_header_norms(FILE *out)
{
	fprintf(out,"# <Time> <NormQP L=0> <NormPhonons L=0> ... <NormQP L=Lmax> <NormPhonons L=Lmax> <TotalNorm> <Re(AlignmentCosine)> <Im(AlignmentCosine)> <Re(Cos^2)> <Im(Cos^2)> <Re(S(t))> <Im(S(t))> <LaserIntensity> <BathIntensity> <RotationalEnergy> <BosonsRotationalEnergy>\n");
}

void print_header_summary(FILE *out)
{
	fprintf(out,"# <Time> <LaserIntensity> <BathIntensity> <Norm> <NormQP> <NormPhonons> <Re(AlignmentCosine)> <Im(AlignmentCosine)> <Re(S(t))> <Im(S(t))> <Completion%%> <LocalNormErr> <TimeDE> <TimeCos> <Torque> <RotationalEnergy> <RotationalEnergy> <BosonsRotationalEnergy>\n");
}

void dump_phonons(FILE *out,struct bigpsi_t *psi,int L,struct configuration_t *config)
{
	int offset=L*(2+10*config->gridpoints);
	double *y=&psi->y[offset];

	for(int d=0;d<config->gridpoints;d++)
	{
        	double gridstep=config->cutoff/config->gridpoints;
		double k=d*gridstep;
		double ti=psi->t;

		double complex phase2m2,phase2m1,phase20,phase21,phase22;
		double complex alpha2m2,alpha2m1,alpha20,alpha21,alpha22;

		phase2m2=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),ti,config);
		phase2m1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),ti,config);
		phase20=timephase(-(L*(L+1.0f)+omegak(k,config)+6.0f),ti,config);
		phase21=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),ti,config);
		phase22=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),ti,config);

		alpha2m2=phase2m2*(y[2+10*d]+I*y[2+10*d+1]);
		alpha2m1=phase2m1*(y[2+10*d+2]+I*y[2+10*d+3]);
		alpha20=phase20*(y[2+10*d+4]+I*y[2+10*d+5]);
		alpha21=phase21*(y[2+10*d+6]+I*y[2+10*d+7]);
		alpha22=phase22*(y[2+10*d+8]+I*y[2+10*d+9]);

		fprintf(out,"%f %f %f %f %f %f %f %f %f %f %f %f\n",ti,k,creal(alpha2m2),cimag(alpha2m2),creal(alpha2m1),cimag(alpha2m1),creal(alpha20),cimag(alpha20),creal(alpha21),cimag(alpha21),creal(alpha22),cimag(alpha22));
	}
}

int do_single(struct configuration_t *config)
{
	struct bigpsi_t *psi;

	FILE **outgs;
	FILE **outalphas;
	FILE *norms;

	char fname[1024];
	int timedivs;

	double *y0,t0;
	bool overlap_snapshot_saved=false;

	outgs=malloc(sizeof(FILE *)*config->maxl);
	outalphas=malloc(sizeof(FILE *)*config->maxl);

	psi=bigpsi_init(config,BIGPSI_INIT_FROM_CONFIG,0,0);

	for(int c=0;c<config->maxl;c++)
		outgs[c]=outalphas[c]=NULL;
	
	norms=NULL;

	for(int c=0;c<config->maxl;c++)
	{
		snprintf(fname,1024,"%s.outg.%d.dat",config->prefix,c);
		fname[1023]='\0';

		if(!(outgs[c]=fopen_mkdir(fname,"w+")))
		{
			fprintf(stderr,"Couldn't open output file!\n");
			goto cleanup;
		}

		print_header_g(outgs[c],psi->params[c].L,config);

		if(config->writephonons==true)
		{
			snprintf(fname,1024,"%s.outalpha.%d.dat",config->prefix,c);
			fname[1023]='\0';

			if(!(outalphas[c]=fopen_mkdir(fname,"w+")))
			{
				fprintf(stderr,"Couldn't open output file!\n");
				goto cleanup;
			}

			print_header_alpha(outalphas[c],psi->params[c].L,config);
		}
	}

	snprintf(fname,1024,"%s.norms.dat",config->prefix);
	if(!(norms=fopen_mkdir(fname,"w+")))
	{
		fprintf(stderr,"Couldn't open output file!\n");
		goto cleanup;
	}

	print_header_norms(norms);
	print_header_summary(stdout);

	timedivs=(config->endtime-config->starttime)/config->timestep;
	for(int c=0;c<=timedivs;c++)
	{
		double ti=config->starttime+c*(config->endtime-config->starttime)/((double)(timedivs));
		double bath_intensity;
		double completion,previousnorm,elapsed_time1,elapsed_time2;
		double complex ac,cs;
		double complex S=0.0f;
		double trq;
		double rotational_e,bosons_rotational_e;

		struct timeval starttime,endtime;

		/*
			If we are past the the overlap zero-time, and we don't have
			a snapshot it's time to take it....
		*/

		if((overlap_snapshot_saved==false)&&(config->overlap==true))
		{
			if(config->overlapt0<ti)
			{
				int size=(2+10*config->gridpoints)*psi->nrpsis*sizeof(double);
				
				y0=malloc(size);
				memcpy(y0,psi->y,size);
				t0=psi->t;

				overlap_snapshot_saved=true;
			}
		}

		/*
			...on the other hand, if we have a snapshot, we have to calcolate the overlap!
		*/

		if((overlap_snapshot_saved==true)&&(config->overlap==true))
			S=overlapS(psi,y0,t0,config);

		/*
			Real time evolution: at first we calculate the data at each step
		*/

		gettimeofday(&starttime,NULL);
		bigpsi_apply_step(psi,ti,&previousnorm,config);
		gettimeofday(&endtime,NULL);

		elapsed_time1=(endtime.tv_sec-starttime.tv_sec)*1000.0;
		elapsed_time1+=(endtime.tv_usec-starttime.tv_usec)/1000.0;
		elapsed_time1/=1000;

		if(config->ramp)
			bath_intensity=adiabatic_ramp(ti,config);
		else
			bath_intensity=1.0f;

		gettimeofday(&starttime,NULL);
		ac=(config->cos2d==true)?(costheta2d(psi,config)):(0.0f);
		cs=costhetasquared(psi,config);
		trq=torque(psi,config->startl[0],config->startm,config);
		rotational_e=rotational_energy(psi,config);
		bosons_rotational_e=bosons_rotational_energy(psi,config);
		gettimeofday(&endtime,NULL);

		elapsed_time2=(endtime.tv_sec-starttime.tv_sec)*1000.0;
		elapsed_time2+=(endtime.tv_usec-starttime.tv_usec)/1000.0;
		elapsed_time2/=1000;

		/*
			...then we write the output to outgs and outalphas, if requested...
		*/

		for(int n=0;n<config->maxl;n++)
		{
			int L=psi->params[n].L;
			int offset=n*(2+10*config->gridpoints);
			
			double *y=&psi->y[offset];

			fprintf(outgs[n],"%f %f %f %f\n",ti,y[0],y[1],norm_qp(ti,y,&psi->params[n],config));

			if(config->writephonons==true)
			{
				for(int d=0;d<config->gridpoints;d++)
				{
		                	double gridstep=config->cutoff/config->gridpoints;
					double k=d*gridstep;

					double complex phase2m2,phase2m1,phase20,phase21,phase22;
					double complex alpha2m2,alpha2m1,alpha20,alpha21,alpha22;

					//if((d%5)!=0)
					//	continue;

					phase2m2=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),ti,config);
					phase2m1=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),ti,config);
					phase20=timephase(-(L*(L+1.0f)+omegak(k,config)+6.0f),ti,config);
					phase21=timephase(-(L*(L+1.0f)+omegak(k,config)+4.0f),ti,config);
					phase22=timephase(-(L*(L+1.0f)+omegak(k,config)-2.0f),ti,config);

					alpha2m2=phase2m2*(y[2+10*d]+I*y[2+10*d+1]);
					alpha2m1=phase2m1*(y[2+10*d+2]+I*y[2+10*d+3]);
					alpha20=phase20*(y[2+10*d+4]+I*y[2+10*d+5]);
					alpha21=phase21*(y[2+10*d+6]+I*y[2+10*d+7]);
					alpha22=phase22*(y[2+10*d+8]+I*y[2+10*d+9]);

					fprintf(outalphas[n],"%f %f %f %f %f %f %f %f %f %f %f %f\n",ti,k,creal(alpha2m2),cimag(alpha2m2),creal(alpha2m1),cimag(alpha2m1),creal(alpha20),cimag(alpha20),creal(alpha21),cimag(alpha21),creal(alpha22),cimag(alpha22));
				}
				
				/*
					This last newline is needed in order for gnuplot to easily read the output file.
				*/

				fprintf(outalphas[n],"\n");
			}
		}

		for(int n=0;n<config->maxl;n++)
		{
			fflush(outgs[n]);

			if(config->writephonons==true)
				fflush(outalphas[n]);
		}

		/*
			...then we write the output to the 'norms' file...
		*/

		fprintf(norms,"%f ",psi->t);	

		for(int d=0;d<psi->nrpsis;d++)
		{
			int offset=d*(2+10*config->gridpoints);

			fprintf(norms,"%f %f ",norm_qp(ti,&psi->y[offset],&psi->params[d],config),norm_phonons(ti,&psi->y[offset],&psi->params[d],config));
		}

		fprintf(norms,"%f %f %f %f %f %f %f %f %f %f %f\n",total_norm(psi),creal(ac),cimag(ac),creal(cs),cimag(cs),creal(S),cimag(S),get_laser_intensity(config->fluence,config->duration,ti,config),bath_intensity,rotational_e,bosons_rotational_e);
		fflush(norms);

		/*
			...and finally a summary on the standard output.
		*/

		completion=100.0f*(ti-config->starttime)/(config->endtime-config->starttime);

		printf("%f %f %f %f %f %f %f %f %f %f %f %f %f%% ",ti,get_laser_intensity(config->fluence,config->duration,ti,config),bath_intensity,total_norm(psi),total_norm_qp(psi),total_norm_phonons(psi),creal(ac),cimag(ac),creal(cs),cimag(cs),creal(S),cimag(S),completion);
		printf("%f %f %f %f %f %f\n",previousnorm,elapsed_time1,elapsed_time2,trq,rotational_e,bosons_rotational_e);
		fflush(stdout);
	}

	if(config->savefinalconf==true)
	{
		FILE *finalconf;
		
		snprintf(fname,1024,"%s.final.dat",config->prefix);
		if(!(finalconf=fopen_mkdir(fname,"w+")))
		{
			fprintf(stderr,"Couldn't open output file '%s'!\n",fname);
			goto cleanup;
		}

		bigpsi_serialize(psi,finalconf);

		if(finalconf!=NULL)
			fclose(finalconf);
	}

	cleanup:

	bigpsi_fini(psi);

	for(int c=0;c<config->maxl;c++)
	{
		if(outgs[c])
			fclose(outgs[c]);
		
		if(config->writephonons==true)
		{
			if(outalphas[c])
				fclose(outalphas[c]);
		}
	}

	if(norms)
		fclose(norms);

	if(outgs)
		free(outgs);

	if(config->writephonons==true)
		if(outalphas)
			free(outalphas);

	return 0;
}

int do_ini_file(char *inifile)
{
	struct configuration_t config;

	load_config_defaults(&config);
	load_molecules_files("molecules.ini");

	if(ini_parse(inifile,configuration_handler,&config)<0)
	{
		printf("Can't load '%s'\n",inifile);
		return 1;
	}

	save_ini_backup(&config,inifile);

	printf("Angulon dynamics, strong coupling.\n");
	printf("Loaded configuration from: %s\n",inifile);

	printf("\nOutput-related options:\n");
	printf("\tOutput prefix: %s\n",config.prefix);
	printf("\tOutput mode: %s\n",(config.writephonons==false)?("light (no phonons)"):("full (incl. phonons)"));
	printf("\tAlignment cosine: %s\n",(config.cos2d==false)?("not calculated (only 3D cosine)"):("calculated"));
	printf("\tUsing molecular cosine: %s\n",(config.altcos==false)?("false"):("true"));

	if(config.savefinalconf==true)
		printf("\tFinal configuration will be saved to: %s.final.dat\n",config.prefix);

	printf("\nGeneral options:\n");
	printf("\tMaximum L: %d\n",config.maxl);
	printf("\tStart time: %f\n",config.starttime);
	printf("\tEnd time: %f\n",config.endtime);
	printf("\tTimestep: %f\n",config.timestep);
	printf("\tFree evolution: %s\n",(config.freeevolution==true)?("true"):("false"));
	printf("\tBosons at finite temperature: %s\n",(config.bosonsfinitetemperature==true)?("true"):("false"));

	if(config.freeevolution==true)
	{
		/*
			In case of a free evolution we just need a minimal number of points in the grid.
			TODO: Would 0 work?
		*/
		
		config.gridpoints=5;
	}

	printf("\nMolecule type: ");
	
	switch(config.moleculetype)
	{
		case MOLECULE_I2:
		printf("I2\n");
		break;

		case MOLECULE_CS2:
		printf("CS2\n");
		break;

		case MOLECULE_OCS:
		printf("OCS\n");
		break;

		default:
		printf("ERROR!\n");
		break;
	}

	printf("\nInitial conditions:\n");

	if(config.mixture==false)
	{
		printf("\tSingle initial state:\n");
		printf("\tL: ");

		for(int d=0;d<config.nrl;d++)
		{
			printf("%d",config.startl[d]);

			if((d+1)!=config.nrl)
				printf(", ");
		}

		if(config.nrl<=0)
		{
			printf("Error: If not working with a mixture, please specify at least on L.");
			exit(0);
		}

		printf("\n\tM: %d\n",config.startm);
	}
	else
	{
		printf("\tUsing a statistical mixture.\n");		
	}

	if(config.freeevolution==false)
	{
		printf("\nGrid:\n");
		printf("\tk-cutoff: %f\n",config.cutoff);
		printf("\tPoints in the k-grid: %d\n",config.gridpoints);
	}

	printf("\nPotential:\n");
	printf("\tu0: %f\n",config.u0);
	printf("\tu2: %f\n",config.u2);
	printf("\tr0: %f\n",config.r0);
	printf("\tr2: %f\n",config.r2);
	printf("\tDensity: %f\n",config.density);

	if(config.morse==true)
		printf("\tUsing effective Morse potential.\n");
	else
		printf("\tUsing Gaussian form factors.\n");

	printf("\nNumerical precision:\n");
	printf("\thstart: %e\n",config.hstart);
	printf("\tepsabs: %e\n",config.epsabs);
	printf("\tepsrel: %e\n",config.epsrel);
	printf("\tNormalization after each step: %s\n",(config.normalize==true)?("ON"):("OFF"));

	printf("\nDispersion: ");

	if(config.dispersion_is_experimental==true)
		printf("from experimental data.\n");
	else
		printf("linear, with sound speed %f (in units of B).\n",config.soundspeed);

	printf("\nTransformation:\n");

	if(config.wtype==1)
		printf("\tW = \\omega_k\n");
	else if(config.wtype==2)
		printf("\tW = \\omega_k + \\lambda(\\lambda + 1)\n");
	else
		printf("\tError: invalid W transformation specified!\n");

	if(config.fscale==true)
		printf("\tAdditional 'f' scaling: ON\n");
	else
		printf("\tAdditional 'f' scaling: OFF\n");

	printf("\nAdiabatic ramp:\n");
	printf("\tRamp: %s\n",(config.ramp==true)?("ON"):("OFF"));
	printf("\tRamp center: %f\n",config.rampcenter);
	printf("\tRamp delta: %f\n",config.rampdelta);

	printf("\nLaser pulse:\n");
	printf("\tLaser: %s\n",(config.laser==true)?("ON"):("OFF"));
	printf("\tFluence (J/cm^2): %f\n",config.fluence);

	if(config.shapefile!=NULL)
	{
		printf("\tLoading pulse shape from '%s'...",config.shapefile);
		load_almost_gaussian(config.shapefile,&config.shapemax);
		printf(" Done!\n");
	}
	else
	{
		printf("\tDuration (ps, FWHM): %f\n",config.duration);
	}

	printf("\nDynamical overlap S(t) = <psi(t=0)|psi(t)> :\n");
	printf("\tOverlap: %s\n",(config.overlap==true)?("ON"):("OFF"));
	printf("\tReference time: %f\n",config.overlapt0);

	printf("\nStarting simulation...\n");

	if(config.mixture==false)
		do_single(&config);
	else
		do_mixture(&config);

	printf("\nCleaning up... ");

	if(config.shapefile!=NULL)
		unload_almost_gaussian();

	printf("Bye!\n");

	return 0;
}

int main(int argc,char **argv)
{
	if(argc<2)
	{
		printf("Usage: %s <inifile> [<otherinifiles> ...]\n",argv[0]);
		return 0;
	}
	
	for(int c=1;c<argc;c++)
		do_ini_file(argv[c]);
	
	return 0;
}
