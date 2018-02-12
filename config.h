#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

struct configuration_t
{
	char *prefix;
	bool writephonons;
	bool cos2d;
	bool savefinalconf;
	
	int maxl;
	double starttime;
	double endtime;
	double timestep;
	bool freeevolution;

#define MOLECULE_I2	(20)
#define MOLECULE_CS2	(21)

	short moleculetype;

	bool mixture;
	int startl;
	int startm;

	double cutoff;
	int gridpoints;

	double u0;
	double u2;
	double r0;
	double r2;
	double density;

	bool ramp;
	double rampcenter;
	double rampdelta;

	double hstart;
	double epsabs;
	double epsrel;
	bool normalize;

	bool laser;
	double duration;
	double milliwatts;
	char *shapefile;
};

int configuration_handler(void *user,const char *section,const char *name,const char *value);

#endif //__CONFIG_H__