#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

struct configuration_t
{
	char *prefix;
	bool writephonons;
	bool cos2d;
	bool savefinalconf;
	bool altcos;

	int maxl;
	double starttime;
	double endtime;
	double timestep;
	bool freeevolution;

	short wtype;
	bool fscale;

#define MOLECULE_I2	(20)
#define MOLECULE_CS2	(21)

	short moleculetype;

	bool mixture;
	int nrl;
	int startl[128];
	int startm;

	double cutoff;
	int gridpoints;

	double u0;
	double u2;
	double r0;
	double r2;
	double density;
	bool morse;

	bool ramp;
	double rampcenter;
	double rampdelta;

	double hstart;
	double epsabs;
	double epsrel;
	bool normalize;

	bool laser;
	double duration;
	double fluence;
	char *shapefile;
	
	bool overlap;
	double overlapt0;
};

int configuration_handler(void *user,const char *section,const char *name,const char *value);

#endif //__CONFIG_H__
