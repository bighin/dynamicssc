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

void load_almost_gaussian(char *fname,double *max)
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
	*max=0.0f;

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
		
		if(b>*max)
			*max=b;
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

double get_laser_intensity(double fluence,double pulse_duration,double t,struct configuration_t *config)
{
	double etamax,pulse_fwhm,pulse_sigma;

	/*
		The FWHM is the quantity reported in experiments. At first we convert it in
		units of B, so that:

		- pulse_duration is the FWHM in ps, as reported in experiments
		- pulse_fwhm is the FWHM in units of B

		and finally pulse_sigma is the usually \sigma of the Gaussian shape, in units of B.
	*/

	pulse_fwhm=pulse_duration/B_in_ps(config);
	pulse_sigma=pulse_fwhm/2.35;

	/*
		The formula we use is:
	
		eta(t) = \eta_max sqrt(2 pi \sigma^2) f(t)
	
		where f(t) is a normalized Gaussian and
	
		\eta_max = 9.91175 * \Delta [Å^3] F [J/cm^2] / ( B [cm^-1] FWHM [ps])
	
		Note: following Bretislav's book one would get 9.35 as conversion factor,
		a good approximation, but not accurate enough in this case.
	*/

	switch(config->moleculetype)
	{
		case MOLECULE_I2:
		etamax=9.91175*(6.0993/0.03739/pulse_duration)*fluence;
		break;

		case MOLECULE_CS2:
		etamax=9.91175*(10.3/0.10901/pulse_duration)*fluence;
		break;

		case MOLECULE_OCS:
		etamax=9.91175*(4.67/0.20286/pulse_duration)*fluence;
		break;

		default:
		fprintf(stderr,"Fatal error: unknown molecular species!\n");
		exit(0);
	}

	etamax*=pulse_sigma*sqrt(2*M_PI);

	if(fabs(t)>=5.0f*pulse_sigma)
		return 0.0f;

	/*
		This is tricky!
	
		The following assumes that the loaded pulse shape is normalized to
		unity, so that the integral of \eta(t) of time is the same as using
		a Gaussian shape.
	
		This is the right thing to do, since the user is specifying the fluence.
	
		However, this means that etamax is the maximum eta of an equivalent
		Gaussian-shaped pulse, but is no longer the *actual* maximum eta.
	*/

	if(config->shapefile!=NULL)
		return etamax*almost_gaussian(t,config);

	return etamax*gaussian(t,0.0f,pulse_sigma);
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

		case MOLECULE_OCS:
	        return 26.1;
		break;
	}

	fprintf(stderr,"Fatal error: unknown molecular species!\n");
	exit(0);
}
