#ifndef __CONFIG_H__
#define __CONFIG_H__

#include <stdbool.h>

#ifndef GITCOMMIT
#define GITCOMMIT "unknown"
#endif

struct configuration_t
{
	char *prefix;
	char *inipath;
	bool writephonons;
	bool cos2d;
	bool savefinalconf;
	bool altcos;

	int maxl;
	double starttime;
	double endtime;
	double timestep;

#define EVOLUTION_FREE		(39)
#define EVOLUTION_1PHONON	(40)
#define EVOLUTION_1PHONONFT	(41)
#define EVOLUTION_COHERENT	(42)

	int evolution;
	double temperature;

	bool dispersion_is_experimental;
	double soundspeed;

	short wtype;
	bool fscale;

#define MOLECULE_I2	(20)
#define MOLECULE_CS2	(21)
#define MOLECULE_OCS	(22)

	short moleculetype;
        bool centrifugal;
        double centrifugalD;
	int centrifugalLcutoff;

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
	double shapemax;
	
	bool overlap;
	double overlapt0;
};

void load_config_defaults(struct configuration_t *config);
int configuration_handler(void *user,const char *section,const char *name,const char *value);
void save_ini_backup(struct configuration_t *config,char *inifile);

#endif //__CONFIG_H__
