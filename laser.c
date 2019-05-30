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
double almost_gaussian_min_x,almost_gaussian_max_x,almost_gaussian_peak;

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

	almost_gaussian_peak=almost_gaussian_y[0];
	for(int c=1;c<cnt;c++)
		if(almost_gaussian_y[c]>almost_gaussian_peak)
			almost_gaussian_peak=almost_gaussian_y[c];

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
	/*
		The innermost conversion factor is needed for converting between
		units of B and ps, while the outermost one is needed to keep the
		normalisation of overall curve.

		Indeed the idea of the following function get_laser_intensity()
		is that almost_gaussian() (and not almost_gaussian_ps()!) looks
		pretty much like a normalised Gaussian.
	*/

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

	etamax=9.91175*(config->Delta_alpha/config->B_in_cms_minus_one/pulse_duration)*fluence;

	etamax*=pulse_sigma*sqrt(2*M_PI);

	if(fabs(t)>=5.0f*pulse_sigma)
		return 0.0f;

	/*
		This is tricky!
	
		The following assumes that the loaded pulse shape is normalized to
		unity, so that the integral of \eta(t) with respect to time is the same
		as using a Gaussian shape.
	
		This is the right thing to do, since the user is specifying the fluence.

		However, this means that etamax is the maximum eta of an equivalent
		Gaussian-shaped pulse, but is no longer the *actual* maximum eta.
	*/

	if(config->shapefile!=NULL)
		return etamax*almost_gaussian(t,config);

	return etamax*gaussian(t,0.0f,pulse_sigma);
}

double get_peak_intensity(double fluence,double pulse_duration,struct configuration_t *config)
{
	double etamax,pulse_fwhm,pulse_sigma;

	/*
		Same reasoning as in the function above, but here
		we calculate the *peak* intensity.
	*/

	pulse_fwhm=pulse_duration/B_in_ps(config);
	pulse_sigma=pulse_fwhm/2.35;

	etamax=9.91175*(config->Delta_alpha/config->B_in_cms_minus_one/pulse_duration)*fluence;
	etamax*=pulse_sigma*sqrt(2*M_PI);

	/*
		The need for a B_in_ps() factor here is clear looking
		at the comment in almost_gaussian().
	*/

	if(config->shapefile!=NULL)
		return etamax*B_in_ps(config)*almost_gaussian_peak;

	return etamax*gaussian(0.0,0.0f,pulse_sigma);
}

double B_in_ps(struct configuration_t *config)
{
	double B_in_THz=config->B*0.001f;

	/*
		Mind the 2 \pi factor!
	*/

	return 1.0f/(2.0f*M_PI*B_in_THz);
}
