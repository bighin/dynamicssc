#ifndef __LASER_H__
#define __LASER_H__

#include "auxx.h"

double gaussian(double x,double mu,double sigma);

struct interpolation_t *almost_gaussian_interpolation;
double almost_gaussian_min_x,almost_gaussian_max_x;

void load_almost_gaussian(char *fname,double *max);
void unload_almost_gaussian(void);
double almost_gaussian_ps(double x);
double almost_gaussian(double x,struct configuration_t *config);

double get_laser_intensity(double mw,double pulse_duration,double t,struct configuration_t *config);

double B_in_ps(struct configuration_t *config);

#endif //__LASER_H__
