#include <stdlib.h>
#include <string.h>

#include "molecules.h"
#include "inih/ini.h"

int find_molecule_id(struct molecule_db_t *moldb,const char *section)
{
	for(int c=0;c<moldb->index;c++)
		if(strcmp(moldb->ids[c],section)==0)
			return c;

	if((moldb->index+1)>=MAX_NR_OF_MOLECULES)
	{
		printf("Error: maximum number of molecules in database reached.\n");
		exit(0);
	}

	moldb->ids[moldb->index]=malloc(sizeof(char)*128);
	snprintf(moldb->ids[moldb->index],128,"%s",section);
	moldb->ids[moldb->index][127]='\0';

	return ++moldb->index;
}

int molecules_configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct molecule_db_t *moldb=(struct molecule_db_t *)(user);

	int id=find_molecule_id(moldb,section);

	if(strcmp(name,"A")==0)
		moldb->As[id]=atof(value);

	else if(strcmp(name,"B")==0)
		moldb->Bs[id]=atof(value);

	else if(strcmp(name,"alpha_par_volume")==0)
		moldb->alphapars[id]=atof(value);

	else if(strcmp(name,"alpha_perp_volume")==0)
		moldb->alphaperps[id]=atof(value);

	else if(strcmp(name,"odd")==0)
		moldb->odds[id]=atoi(value);

	else if(strcmp(name,"even")==0)
		moldb->evens[id]=atoi(value);

	return 1;
}

void load_molecules_files(char *molfile)
{
	struct molecule_db_t moldb;

	moldb.index=0;
	
	if(ini_parse(molfile,molecules_configuration_handler,&moldb)<0)
	{
		printf("Can't load '%s'\n",molfile);
		return;
	}
	
	for(int c=0;c<moldb.index;c++)
	{
		printf("Molecule: %s\n",moldb.ids[c]);

		printf("\tA: %f\n",moldb.As[c]);
		printf("\tB: %f\n",moldb.Bs[c]);
		printf("\talpha_par: %f\n",moldb.alphapars[c]);
		printf("\talpha_perp: %f\n",moldb.alphaperps[c]);
		printf("\tEven abundance: %d\n",moldb.evens[c]);
		printf("\tOdd abundance: %d\n",moldb.odds[c]);
	
		if((c+1)!=moldb.index)
			printf("\n");
	}
}
