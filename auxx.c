#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include "auxx.h"
#include "spline.h"

/*
	Spline interpolation, given a grid of points, essentially a nice wrapper
	over the code contained in Numerical Recipes (NR).

	Note that the code in NR is __not__ thread safe, the static keyword has to
	be removed in the splint() function.
*/

struct interpolation_t *init_interpolation(double *x,double *y,int n)
{
	struct interpolation_t *ret;
	
	if(!(ret=malloc(sizeof(struct interpolation_t))))
		return NULL;

	if(!(ret->x=malloc(sizeof(double)*n)))
	{
		if(ret)
			free(ret);
	
		return NULL;
	}

	if(!(ret->y=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			free(ret);
		}
		
		return NULL;
	}

	if(!(ret->y2=malloc(sizeof(double)*n)))
	{
		if(ret)
		{
			if(ret->x)
				free(ret->x);

			if(ret->y)
				free(ret->y);
		
			free(ret);
		}

		return NULL;
	}

	memcpy(ret->x,x,sizeof(double)*n);
	memcpy(ret->y,y,sizeof(double)*n);
	ret->n=n;
	
	spline(ret->x,ret->y,ret->n,0.0f,0.0f,ret->y2);

	return ret;
}

void fini_interpolation(struct interpolation_t *it)
{
	if(it)
	{
		if(it->x)
			free(it->x);

		if(it->y)
			free(it->y);

		if(it->y2)
			free(it->y2);
	
		free(it);
	}
}

void copy_interpolation(struct interpolation_t *dst,struct interpolation_t *src)
{
	assert(dst);
	assert(src);

	dst->n=src->n;

	memcpy(dst->x,src->x,sizeof(double)*dst->n);
	memcpy(dst->y,src->y,sizeof(double)*dst->n);
	memcpy(dst->y2,src->y2,sizeof(double)*dst->n);
}

double get_point(struct interpolation_t *it,double x)
{
	double y;
	
	splint(it->x,it->y,it->y2,it->n,x,&y);
	
	return y;
}

FILE *fopen_mkdir(const char *name,const char *mode)
{
	char *mname=strdup(name);

	for(int i=0;mname[i]!='\0';i++)
	{
		if((i>0)&&((mname[i]=='\\')||(mname[i]=='/')))
		{
			char slash=mname[i];

			mname[i]='\0';

			if(!mkdir(mname,S_IRWXU|S_IRWXG|S_IRWXO))
			{
				free(mname);
				return NULL;
			}
			
			mname[i]=slash;
		}
	}

	free(mname);
	return fopen(name,mode);
}

double arccot(double x)
{
	return (M_PI/2.0f)-atan(x);
}
