#ifndef __AUX_H__
#define __AUX_H__

#include <stdio.h>
#include <math.h>

struct interpolation_t
{
	double *x;
	double *y;

	int n;
	
	double *y2;
};

struct interpolation_t *init_interpolation(double *x,double *y,int n);
void fini_interpolation(struct interpolation_t *it);
void copy_interpolation(struct interpolation_t *dst,struct interpolation_t *src);
double get_point(struct interpolation_t *it,double x);

void fcopy(FILE *f1,FILE *f2);
FILE *fopen_mkdir(const char *name,const char *mode);

double arccot(double x);

#endif //__AUX_H__
