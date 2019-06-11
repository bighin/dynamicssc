#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <stdbool.h>

#include "inih/ini.h"
#include "config.h"
#include "auxx.h"
#include "molecules.h"

void load_config_defaults(struct configuration_t *config)
{
	memset(config,0,sizeof(struct configuration_t));

	config->prefix="dsc";
	config->writephonons=false;
	config->cos2d=false;
	config->savefinalconf=false;
	config->altcos=false;

	config->maxl=32;
	config->starttime=-5.0f;
	config->endtime=5.0f;
	config->timestep=0.01f;
	config->evolution=EVOLUTION_FREE;
	config->temperature=0.38;

	config->dispersion_is_experimental=true;
	config->soundspeed=23.00;

	config->wtype=1;
	config->fscale=true;

	config->moleculetype=0;
        config->centrifugal=false;
        config->centrifugalD=0.0f;
        config->centrifugalLcutoff=20;
	config->highercorrection=false;
        config->realspectrum=false;

	config->mixture=false;
	config->nrl=1;
	config->startl[0]=0;
	config->startm=0;

	config->cutoff=8.0f;
	config->gridpoints=800;

	config->u0=0.0f;
	config->u2=600.0f;
	config->r0=0.0f;
	config->r2=1.0f;
	config->density=74;
	config->morse=false;

	config->ramp=false;
	config->rampcenter=0.0f;
	config->rampdelta=2.0f;

	config->hstart=1e-6;
	config->epsabs=1e-6;
	config->epsrel=1e-6;
	config->normalize=true;

	config->laser=true;
	config->duration=20;
	config->fluence=0.4;
	config->shapefile=NULL;
	config->shapemax=1.0f;

	config->overlap=false;
	config->overlapt0=-2.5;

	config->mixture_nr_states=0;
}

int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);
	struct molecule_db_t *moldb=pconfig->moldb;

#define MATCH(s,n) (strcmp(section,s)==0)&&(strcmp(name,n)==0)

	if(MATCH("output","prefix"))
	{
		if(strcmp(value,"auto")==0)
		{
			if(!strstr(pconfig->inipath,".ini"))
			{
				printf("Error: using automatic prefix, but the configuration file path does not contain '.ini'\n");
				exit(0);
			}

			pconfig->prefix=find_and_replace(pconfig->inipath,".ini","");
		}
		else
		{
			pconfig->prefix=strdup(value);
		}
	}
	if(MATCH("output","writephonons"))
	{
		if(!strcasecmp(value,"true"))
			pconfig->writephonons=true;
		else
			pconfig->writephonons=false;
	}
	if(MATCH("output","cos2d"))
	{
		if(!strcasecmp(value,"true"))
			pconfig->cos2d=true;
		else
			pconfig->cos2d=false;
	}
	if(MATCH("output","savefinalconf"))
	{
		if(!strcasecmp(value,"true"))
			pconfig->savefinalconf=true;
		else
			pconfig->savefinalconf=false;
	}
	if(MATCH("output","altcos"))
	{
		if(!strcasecmp(value,"true"))
			pconfig->altcos=true;
		else
			pconfig->altcos=false;
	}
	else if(MATCH("general","maxl"))
	{
		pconfig->maxl=atoi(value);
	}
	else if(MATCH("general","starttime"))
	{
		pconfig->starttime=atof(value);
	}
	else if(MATCH("general","endtime"))
	{
		pconfig->endtime=atof(value);
	}
	else if(MATCH("general","timestep"))
	{
		pconfig->timestep=atof(value);
	}
	else if(MATCH("general","freeevolution"))
	{
		fprintf(stderr,"Error: you are using the old option 'freeevolution', please use the new one\n");
		fprintf(stderr,"'evolution' which can take the values: 'free', '1phonon', '1phononft' or 'coherent'.\n");
		exit(0);
	}
	else if(MATCH("general","evolution"))
	{
		if(!strcasecmp(value,"free"))
		{
			pconfig->evolution=EVOLUTION_FREE;
		}
		else if(!strcasecmp(value,"1phononft"))
		{
			pconfig->evolution=EVOLUTION_1PHONONFT;
		}
		else if(!strcasecmp(value,"1phonon"))
		{
			pconfig->evolution=EVOLUTION_1PHONON;
		}
		else if(!strcasecmp(value,"coherent"))
		{
			pconfig->evolution=EVOLUTION_COHERENT;
		}
		else
		{
			fprintf(stderr,"Error: unknown evolution type (%s)\n",value);
			exit(0);
		}
	}
	else if(MATCH("general","bosonsfinitetemperature"))
	{
		fprintf(stderr,"Error: you are using the old option 'bosonsfinitetemperature', please use the new one\n");
		fprintf(stderr,"'evolution' which can take the values: 'free', '1phonon', '1phononft' or 'coherent'.\n");
		exit(0);
	}
	else if(MATCH("general","temperature"))
	{
		pconfig->temperature=atof(value);
	}
	else if(MATCH("dispersion","dispersion"))
	{
		if(!strcasecmp(value,"linear"))
			pconfig->dispersion_is_experimental=false;
		else
			pconfig->dispersion_is_experimental=true;
	}
	else if(MATCH("dispersion","soundspeed"))
	{
		pconfig->soundspeed=atof(value);
	}
	else if(MATCH("transformation","wtype"))
	{
		pconfig->wtype=atoi(value);
	}
	else if(MATCH("transformation","fscale"))
	{
		if(!strcasecmp(value,"true"))
			pconfig->fscale=true;
		else
			pconfig->fscale=false;
	}
	else if(MATCH("molecule","type"))
	{
		pconfig->moleculetype=find_molecule_id(moldb,value);

		pconfig->B=moldb->Bs[pconfig->moleculetype];
		pconfig->B_in_cms_minus_one=pconfig->B/29.9792458;
		pconfig->Delta_alpha=moldb->alphapars[pconfig->moleculetype]-moldb->alphaperps[pconfig->moleculetype];

		if(pconfig->moleculetype==-1)
		{
			fprintf(stderr,"Error: invalid molecule type '%s'\n",value);
			exit(0);
		}
	}
	else if(MATCH("molecule","centrifugalD"))
	{
		pconfig->centrifugal=true;
		pconfig->centrifugalD=atof(value);
        }
	else if(MATCH("molecule","centrifugalLcutoff"))
	{
		pconfig->centrifugal=true;
		pconfig->centrifugalLcutoff=atoi(value);
        }
	else if(MATCH("molecule","realspectrum"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->realspectrum=true;
		else
			pconfig->realspectrum=false;
	}
	else if(MATCH("molecule","highercorrection"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->highercorrection=true;
		else
			pconfig->highercorrection=false;
	}
	else if(MATCH("initialconditions","type"))
	{
		if(!strcmp(value,"singlestate"))
			pconfig->mixture=false;
		else if(!strcmp(value,"mixture"))
			pconfig->mixture=true;
		else
			fprintf(stderr,"Warning: invalid initial conditions '%s'\n",value);
	}
	else if(MATCH("initialconditions","l"))
	{
		char *token,*string,*tofree;

		pconfig->nrl=0;		
		tofree=string=strdup(value);

		while((token=strsep(&string,","))!=NULL)
		{
			int lval=atoi(token);

			pconfig->startl[pconfig->nrl]=lval;
			pconfig->nrl++;
		}

		free(tofree);
	}
	else if(MATCH("initialconditions","m"))
	{
		pconfig->startm=atoi(value);
	}
	else if(MATCH("grid","cutoff"))
	{
		pconfig->cutoff=atof(value);
	}
	else if(MATCH("grid","gridpoints"))
	{
		pconfig->gridpoints=atoi(value);
	}
	else if(MATCH("potential","u0"))
	{
		pconfig->u0=atof(value);
	}
	else if(MATCH("potential","u2"))
	{
		pconfig->u2=atof(value);
	}
	else if(MATCH("potential","r0"))
	{
		pconfig->r0=atof(value);
	}
	else if(MATCH("potential","r2"))
	{
		pconfig->r2=atof(value);
	}
	else if(MATCH("potential","density"))
	{
		pconfig->density=atof(value);
	}
	else if(MATCH("potential","morse"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->morse=true;
		else
			pconfig->morse=false;
	}
	else if(MATCH("precision","hstart"))
	{
		pconfig->hstart=atof(value);
	}
	else if(MATCH("precision","epsabs"))
	{
		pconfig->epsabs=atof(value);
	}
	else if(MATCH("precision","epsrel"))
	{
		pconfig->epsrel=atof(value);
	}
	else if(MATCH("precision","normalize"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->normalize=true;
		else
			pconfig->normalize=false;
	}
	else if(MATCH("adiabaticramp","ramp"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->ramp=true;
		else
			pconfig->ramp=false;
	}
	else if(MATCH("adiabaticramp","rampcenter"))
	{
		pconfig->rampcenter=atof(value);
	}
	else if(MATCH("adiabaticramp","rampdelta"))
	{
		pconfig->rampdelta=atof(value);
	}
	else if(MATCH("pulse","laser"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->laser=true;
		else
			pconfig->laser=false;
	}
	else if(MATCH("pulse","duration"))
	{
		pconfig->duration=atof(value);
	}
	if(MATCH("pulse","shapefile"))
	{
		pconfig->shapefile=strdup(value);
	}
	else if(MATCH("pulse","fluence"))
	{
		pconfig->fluence=atof(value);
	}
	else if(MATCH("overlap","overlap"))
	{
		if((!strcasecmp(value,"true"))||(!strcasecmp(value,"on")))
			pconfig->overlap=true;
		else
			pconfig->overlap=false;
	}
	else if(MATCH("overlap","t0"))
	{
		pconfig->overlapt0=atof(value);
	}
	else
	{
		return 0;  /* unknown section/name, error */
	}

	return 1;
}

void save_ini_backup(struct configuration_t *config,char *inifile)
{
	FILE *in,*out;
	char fname[1024];
	
	snprintf(fname,1024,"%s.inibackup.dat",config->prefix);
	fname[1023]='\0';

	if(!(in=fopen(inifile,"r")))
		return;

	if(!(out=fopen(fname,"w+")))
	{
		if(in)
			fclose(in);
	
		return;
	}

	fprintf(out,"; Binary compiled from commit: %s\n",GITCOMMIT);
	fcopy(in,out);

	if(in)
		fclose(in);

	if(out)
		fclose(out);
}
