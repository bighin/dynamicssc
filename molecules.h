#ifndef __MOLECULES_H__
#define __MOLECULES_H__

#include "config.h"

#define MAX_NR_OF_MOLECULES	(128)

struct molecule_db_t
{
	char *ids[MAX_NR_OF_MOLECULES];
	
	double As[MAX_NR_OF_MOLECULES];
	double Bs[MAX_NR_OF_MOLECULES];
	double alphapars[MAX_NR_OF_MOLECULES];
	double alphaperps[MAX_NR_OF_MOLECULES];

	int odds[MAX_NR_OF_MOLECULES];
	int evens[MAX_NR_OF_MOLECULES];

	int index;
};

int find_molecule_id(struct molecule_db_t *moldb,const char *section);
int find_or_add_molecule_id(struct molecule_db_t *moldb,const char *section);
struct molecule_db_t *load_molecules_files(char *molfile,bool verbose);
void fini_moldb(struct molecule_db_t *moldb);

#endif //__MOLECULES_H__
