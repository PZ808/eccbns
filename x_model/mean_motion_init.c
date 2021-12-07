/*
 * mean_motion_init.c
 * $Id: mean_motion_init.c,v 1.2 2009/04/02 13:43:28 pzimmerman Exp $
 */

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_math.h> 
#include "mean_motion_init.h"

double n_init( double n, void *params )
{
	const double t_sun = 4.92549095e-6;

	struct ic_params *p = (struct ic_params *) params;

	double f_init = p->f_init;
	double M   	  = p->M;
	double eta	  = p->eta;
	double e	    = p->e;
	int pn_order  = p->pn_order;
  double n_SI;

  if ( pn_order == 0 ) 
	{
		n_SI = M_PI * f_init - n;
	} 
  else if ( pn_order == 1 )
	{
		n_SI = M_PI * f_init - n * ( 1. + ( 3.* pow(M * n * t_sun , 2./3.) ) / ( 1. - e*e ) ); 
	}
  else if ( pn_order == 2 )
	{
		n_SI = M_PI * f_init - n * ( 1. + ( 3.* pow(M * n * t_sun , 2./3.)) / (1. - e*e ) 
				+ ( pow(M * n * t_sun, 4./3.) / (4.* pow( 1. - e*e, 2.)) ) 
				* (78. - 28.* eta + (51. - 26.* eta) * e*e) );
	}
	else 
	{
		fprintf( stderr, "error in pn order\n" );
		exit( 1 );
	}
  return n_SI;
}  
