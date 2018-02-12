#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include "inih/ini.h"
#include "config.h"

int configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct configuration_t *pconfig=(struct configuration_t *)(user);

#define MATCH(s,n) (strcmp(section,s)==0)&&(strcmp(name,n)==0)

	if(MATCH("output","prefix"))
	{
		pconfig->prefix=strdup(value);
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
		if(!strcasecmp(value,"true"))
			pconfig->freeevolution=true;
		else
			pconfig->freeevolution=false;
	}
	else if(MATCH("transformation","wtype"))
	{
		pconfig->wtype=atoi(value);
	}
	else if(MATCH("molecule","type"))
	{
		if(!strcmp(value,"I2"))
			pconfig->moleculetype=MOLECULE_I2;
		else if(!strcmp(value,"CS2"))
			pconfig->moleculetype=MOLECULE_CS2;
		else
			fprintf(stderr,"Warning: invalid molecule type '%s'\n",value);
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
		pconfig->startl=atoi(value);
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
	else if(MATCH("pulse","milliwatts"))
	{
		pconfig->milliwatts=atof(value);
	}
	else
	{
		return 0;  /* unknown section/name, error */
	}

	return 1;
}
