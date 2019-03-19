#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "molecules.h"
#include "config.h"
#include "inih/ini.h"

int find_or_add_molecule_id(struct molecule_db_t *moldb,const char *section)
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

int find_molecule_id(struct molecule_db_t *moldb,const char *section)
{
	for(int c=0;c<moldb->index;c++)
		if(strcmp(moldb->ids[c],section)==0)
			return c;
	
	return -1;
}

int molecules_configuration_handler(void *user,const char *section,const char *name,const char *value)
{
	struct molecule_db_t *moldb=(struct molecule_db_t *)(user);

	int id=find_or_add_molecule_id(moldb,section);

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

struct molecule_db_t *load_molecules_files(char *molfile,bool verbose)
{
	struct molecule_db_t *moldb=malloc(sizeof(struct molecule_db_t));

	moldb->index=0;
	
	if(ini_parse(molfile,molecules_configuration_handler,moldb)<0)
	{
		printf("Can't load '%s'\n",molfile);
		
		if(moldb)
			free(moldb);
		
		return NULL;
	}

	if(verbose)
	{
		for(int c=0;c<moldb->index;c++)
		{
			printf("Molecule: %s\n",moldb->ids[c]);

			printf("\tA: %f\n",moldb->As[c]);
			printf("\tB: %f\n",moldb->Bs[c]);
			printf("\talpha_par: %f\n",moldb->alphapars[c]);
			printf("\talpha_perp: %f\n",moldb->alphaperps[c]);
			printf("\tEven abundance: %d\n",moldb->evens[c]);
			printf("\tOdd abundance: %d\n",moldb->odds[c]);
		
			if((c+1)!=moldb->index)
				printf("\n");
		}
	}

	return moldb;
}

void fini_moldb(struct molecule_db_t *moldb)
{
	if(moldb)
		free(moldb);
}

void too_many_states(int max_nr_of_states,double threshold)
{
	fprintf(stderr,"The maximum number of states (MIXTURE_MAX_NR_STATES=%d) has been reached before hitting the specified threshold threshold (%f).\n",max_nr_of_states,threshold);
	fprintf(stderr,"Try increasing MIXTURE_MAX_NR_STATES or decreasing the threshold.\n");
}

/*
	Here we use the molecular data to calculate the statistical weights for the mixture
*/

#define MAXJ	(32)

void calculate_mixture_weights(struct molecule_db_t *moldb,struct configuration_t *config,bool verbose)
{
	double weights[MAXJ+1][2*(MAXJ+1)+1],totaleven,totalodd,cumulative,BinK,rcr,threshold;
	int id,maxjeven,maxjodd;

	totaleven=totalodd=0.0f;
	maxjeven=maxjodd=0;

	/*
		We try and find the molecular id in the database
	*/

	id=-1;
	switch(config->moleculetype)
	{
		case MOLECULE_I2:
		id=find_molecule_id(moldb,"I2");
		break;

		case MOLECULE_CS2:
		id=find_molecule_id(moldb,"CS2");
		break;

		case MOLECULE_OCS:
		id=find_molecule_id(moldb,"OCS");
		break;

		default:
		fprintf(stderr,"Fatal error: unknown molecular species!\n");
		exit(0);
	}

	if(id==-1)
	{
		fprintf(stderr,"Fatal: molecules.conf does not contain information for the specified molecule.\n");
		exit(0);		
	}

	/*
		Conversion from GHz to K.
	*/

	BinK=moldb->Bs[id]/20.8366176361328;

	/*
		The rotational constant renormalisation for different molecules.

		Note that the value is hardcoded is (as calculated from the static theory)
		just for calculating the mixture weight. The equations of motion know nothing
		about the rotational constant renormalisation, that appears 'spontaneously'.
	*/

	rcr=1.0f;

	if(config->evolution!=EVOLUTION_FREE)
	{
		switch(config->moleculetype)
		{
			case MOLECULE_I2:
			rcr=0.6f;
			break;

			case MOLECULE_CS2:
			rcr=0.3f;
			break;

			case MOLECULE_OCS:
			rcr=0.36f;
			break;

			default:
			fprintf(stderr,"Fatal error: unknown molecular species!\n");
			exit(0);
		}
	}

	/*
		We have some magic values for the threshold, that end up using a decent number
		of states. This could be improved
	*/

	switch(config->moleculetype)
	{
		case MOLECULE_I2:
		threshold=(config->evolution==EVOLUTION_FREE)?(0.995):(0.975);
		break;

		case MOLECULE_CS2:
		threshold=(config->evolution==EVOLUTION_FREE)?(0.999):(0.999);
		break;

		case MOLECULE_OCS:
		threshold=(config->evolution==EVOLUTION_FREE)?(0.99999):(0.995);
		break;

		default:
		fprintf(stderr,"Fatal error: unknown molecular species!\n");
		exit(0);
	}

	/*
		Finally, we calculate the Boltzmann factors...
	*/

	for(int j=0;j<=MAXJ;j++)
	{
		for(int m=-j;m<=j;m++)
		{
			double weight;
			bool iseven=((j%2)==0)?(true):(false);

			weight=0.0f;

			if((iseven==true)&&(moldb->evens[id]!=0))
				weight=exp(-BinK*rcr*j*(j+1)/config->temperature);

			if((iseven==false)&&(moldb->odds[id]!=0))
				weight=exp(-BinK*rcr*j*(j+1)/config->temperature);

			weights[j][j+m]=weight;

			if(iseven==true)
				totaleven+=weight;
			else
				totalodd+=weight;
			
		}
	}

	/*
		...and we normalise them separatedly for even and odd states.
	*/

	if(moldb->evens[id]!=0)
	{
		cumulative=0.0f;
		for(int j=0;j<=MAXJ;j+=2)
		{
			for(int m=-j;m<=j;m++)
				cumulative+=weights[j][j+m];

			if(cumulative>=threshold*totaleven)
			{
				maxjeven=j;
				break;
			}
		}

		for(int j=0;j<=maxjeven;j+=2)
		{
			for(int m=-j;m<=j;m++)
			{
				weights[j][j+m]/=cumulative;
			}
		}
	}

	if(moldb->odds[id]!=0)
	{
		cumulative=0.0f;
		for(int j=1;j<=MAXJ;j+=2)
		{
			for(int m=-j;m<=j;m++)
				cumulative+=weights[j][j+m];

			if(cumulative>=threshold*totalodd)
			{
				maxjodd=j;
				break;
			}
		}

		for(int j=1;j<=maxjodd;j+=2)
		{
			for(int m=-j;m<=j;m++)
			{
				weights[j][j+m]/=cumulative;
			}
		}
	}

	config->mixture_even_abundance=moldb->evens[id];
	config->mixture_odd_abundance=moldb->odds[id];
	config->mixture_nr_states=0;

#define MAX(a,b)	(((a)>(b))?(a):(b))

	for(int j=0;j<=MAX(maxjeven,maxjodd);j++)
	{
		for(int m=0;m<=j;m++)
		{
			if(m==0)
			{
				if(fabs(weights[j][j+m])>1e-8)
				{
					if(verbose==true)
						printf("State: (j=%d,m=%d) ==> %f\n",j,m,weights[j][j+m]);

					if(config->mixture_nr_states==MIXTURE_MAX_NR_STATES)
						too_many_states(MIXTURE_MAX_NR_STATES,threshold);

					config->mixture_states[config->mixture_nr_states][0]=j;
					config->mixture_states[config->mixture_nr_states][1]=m;
					config->mixture_weights[config->mixture_nr_states]=weights[j][j+m];
					config->mixture_nr_states++;

				}
			}
			else
			{
				if(fabs(weights[j][j+m])>1e-8)
				{
					if(verbose==true)
						printf("State: (j=%d,m=%d) ==> %f\n",j,m,2.0f*weights[j][j+m]);

					if(config->mixture_nr_states==MIXTURE_MAX_NR_STATES)
						too_many_states(MIXTURE_MAX_NR_STATES,threshold);

					config->mixture_states[config->mixture_nr_states][0]=j;
					config->mixture_states[config->mixture_nr_states][1]=m;
					config->mixture_weights[config->mixture_nr_states]=2.0f*weights[j][j+m];
					config->mixture_nr_states++;
				}
			}
		}
	}
}

/*
	Apparently that's what's used in the Aarhus alignment simulator and what we used for
	calculating the OCS weights.

	However this is at odds with the formula in the SM of Jan's paper, which is used to
	implement the previous function.
*/

void calculate_mixture_weights_alt(struct molecule_db_t *moldb,int id,struct configuration_t *config)
{
	double weights[MAXJ+1][2*(MAXJ+1)+1],total,cumulative,BinK,rcr,threshold;
	int maxj;

	threshold=0.999;

	total=0.0f;
	maxj=0;

	/*
		Conversion from GHz to K.
	*/

	BinK=moldb->Bs[id]/20.8366176361328;

	/*
		The rotational constant renormalisation for different molecules.

		Note that the value is hardcoded is (as calculated from the static theory)
		just for calculating the mixture weight. The equations of motion know nothing
		about the rotational constant renormalisation, that appear 'spontaneously'.
	*/

	rcr=1.0f;

	if(config->evolution!=EVOLUTION_FREE)
	{
		switch(config->moleculetype)
		{
			case MOLECULE_I2:
			rcr=0.6f;
			break;

			case MOLECULE_CS2:
			rcr=0.3f;
			break;

			case MOLECULE_OCS:
			rcr=0.36f;
			break;

			default:
			fprintf(stderr,"Fatal error: unknown molecular species!\n");
			exit(0);
		}
	}

	for(int j=0;j<=MAXJ;j++)
	{
		for(int m=-j;m<=j;m++)
		{
			double weight=exp(-BinK*rcr*j*(j+1)/config->temperature);

			weights[j][j+m]=weight;
			total+=weight;
		}
	}

	cumulative=0.0f;
	for(int j=0;j<=MAXJ;j++)
	{
		for(int m=-j;m<=j;m++)
			cumulative+=weights[j][j+m];

		if(cumulative>=threshold*total)
		{
			maxj=j;
			break;
		}
	}

	for(int j=0;j<=maxj;j++)
	{
		for(int m=-j;m<=j;m++)
		{
			weights[j][j+m]/=cumulative;
		}
	}

#define MAX(a,b)	(((a)>(b))?(a):(b))

	for(int j=0;j<=maxj;j++)
	{
		for(int m=0;m<=j;m++)
		{
			if(m==0)
			{
				if(fabs(weights[j][j+m])>1e-8)
					printf("State: (j=%d,m=%d) ==> %f\n",j,m,weights[j][j+m]);
			}
			else
			{
				if(fabs(weights[j][j+m])>1e-8)
					printf("State: (j=%d,m=%d) ==> %f\n",j,m,2.0f*weights[j][j+m]);
			}
		}
	}
}
