/* 
 	eT_3pN_energy.c
 */

#include <stdio.h>
#include <gsl/gsl_math.h>
#include "eT_3pN_energy.h"

double binding_energy_func( double q, void *params )
{
    double t_sun = 4.92549095e-6;

    double *energy_params = (void*) params;

    double M   = energy_params[0];
    double eta = energy_params[1];
	double x;
    double E;
	double xi;
	double z;

	/* 3pN binding energy */
	E = - ( eta * M / 2.) * x * ( 1. - 1./12.* ( 9. + eta ) * x
		- ( 27./8. - 19./8.* eta + 1./24. * eta*eta ) * x*x 
		- ( 675./64. + 35./5184.* eta*eta*eta + 155./96.* eta*eta 
		+ ( 205./96.* M_PI*M_PI - 34445./576.) * eta ) * x*x*x );
	
	z = xi + 2.* E / ( M * eta );
	
	return z;
}
