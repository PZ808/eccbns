/* $Id: eT_ICs.c,v 1.1 2009/02/21 19:02:54 pzimmerman Exp $ */

#include <stdio.h>
#include <gsl/gsl_math.h>

double get_xi_init( double xi, void *params )
{
  double *init_params = (void *)params;
  double M    = init_params[0];
  double eta  = init_params[1];
	double f_gw = init_params[2];
  const double t_sun = 4.92549095e-6;

 	return ( M_PI * f_gw - ( pow(xi, 3./2.) / (M * t_sun) ) * 
 		( 1. + (1./8.) * ( 9. + eta ) * xi
    	+ ( (891./128.) - (201./64.) * eta
    	+	(11./128.) * eta*eta ) * xi*xi
    	+ ( (41445./1024.) - ( (309715./3072.) 
    			- (205./64.) * M_PI*M_PI ) * eta
    			+ (1215./1024.) * eta*eta 
    			+ (45./1024.) * eta*eta*eta ) * xi*xi*xi ) );
}
