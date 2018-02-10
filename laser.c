#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "auxx.h"
#include "config.h"
#include "laser.h"

double gaussian(double x,double mu,double sigma)
{
	return exp(-(x-mu)*(x-mu)/(2.0f*sigma*sigma))/sqrt(2.0f*M_PI*sigma*sigma);
}

struct interpolation_t *almost_gaussian_interpolation;
double almost_gaussian_min_x,almost_gaussian_max_x;

void load_almost_gaussian(char *fname)
{
	FILE *f;
	double *almost_gaussian_x,*almost_gaussian_y;
	int cnt;

	if(!(f=fopen(fname,"r")))
	{
		printf("Error: couldn't open %s for reading!\n",fname);
		exit(0);
	}

	for(cnt=0;(!feof(f))&&(!feof(f));cnt++)
	{
		char line[1024];
		double a,b;
		
		if(!fgets(line,1024,f))
			continue;
		
		if(sscanf(line,"%lf,%lf",&a,&b)!=2)
			continue;
	
		if(cnt==0)
		{
			almost_gaussian_min_x=a;
			almost_gaussian_max_x=a;
		}
	
		if(a<almost_gaussian_min_x)
			almost_gaussian_min_x=a;

		if(a>almost_gaussian_max_x)
			almost_gaussian_max_x=a;
	}

	almost_gaussian_x=malloc(sizeof(double)*cnt);
	almost_gaussian_y=malloc(sizeof(double)*cnt);

	rewind(f);

	for(cnt=0;(!feof(f))&&(!feof(f));cnt++)
	{
		char line[1024];
		double a,b;
		
		if(!fgets(line,1024,f))
			continue;
		
		if(sscanf(line,"%lf,%lf",&a,&b)!=2)
			continue;

		almost_gaussian_x[cnt]=a;
		almost_gaussian_y[cnt]=b;
	}

	almost_gaussian_interpolation=init_interpolation(almost_gaussian_x,almost_gaussian_y,cnt);

	if(almost_gaussian_x)
		free(almost_gaussian_x);

	if(almost_gaussian_y)
		free(almost_gaussian_y);

	if(f)
		fclose(f);
}

void unload_almost_gaussian(void)
{
	fini_interpolation(almost_gaussian_interpolation);
}

double positive_part(double x)
{
	if(x>=0.0f)
		return x;
	
	return 0.0f;
}

double almost_gaussian_ps(double x)
{
	if((x<=almost_gaussian_min_x)||(x>=almost_gaussian_max_x))
		return 0.0f;
	
	return positive_part(get_point(almost_gaussian_interpolation,x));
}

double almost_gaussian(double x,struct configuration_t *config)
{
	return B_in_ps(config)*almost_gaussian_ps(x*B_in_ps(config));
}

double get_laser_intensity(double mw,double pulse_duration,double t,struct configuration_t *config)
{
	double eta,pulse_fwhm,pulse_sigma;
	
	/*
		The fluence is 0.71 for a pulse at 10 mW and increases linearly
		with the pulse intensity.

		The usual formula is valid for I2, to get the same formula for CS2 we
		just multiply it by the polarizability ratio.
	
		This was confirmed by Lars by email.
	*/

	switch(config->moleculetype)
	{
		case MOLECULE_I2:
	        eta=10.5*0.71*(mw/10.0f);
		break;

		case MOLECULE_CS2:
	        eta=(10.5/6.1)*10.5*0.71*(mw/10.0f);
		break;
		
		default:
		fprintf(stderr,"Fatal error: unknown molecular species!\n");
		exit(0);
	}

	/*
		The FWHM is the quantity reported in experiments. At first we convert it in
		units of B, so that:
	
		- pulse_duration is the FWHM in ps, as reported in experiments
		- pulse_fwhm is the FWHM in units of B
	
		and finally pulse_sigma is the usually \sigma of the Gaussian shape, in units of B.
	*/

	pulse_fwhm=pulse_duration/B_in_ps(config);
	pulse_sigma=pulse_fwhm/2.35;

	if(fabs(t)>=5.0f*pulse_sigma)
		return 0.0f;

	if(config->shapefile!=NULL)
		return eta*almost_gaussian(t,config);

	return eta*gaussian(t,0.0f,pulse_sigma);
}

double B_in_ps(struct configuration_t *config)
{
	/*
		This is one of the two numbers that depend on the molecular species --
        	the other being the relation between the fluence and eta -- the factor is:

        	B = 142.0 ps for for I2.
	*/

	switch(config->moleculetype)
	{
		case MOLECULE_I2:
	        return 142.0f;
		break;

		case MOLECULE_CS2:
	        return 48.2;
		break;
	}

	fprintf(stderr,"Fatal error: unknown molecular species!\n");
	exit(0);
}
