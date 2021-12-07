/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *											 *	
 * toy_ode_system.c									 *	
 *											 *
 *    Program sets up the framework for the ode system corresponding to the "toy"  	 *
 *    initial value problem. Flexibility is restricted by the function init_y, which	 *
 *    hardwires the initial time. 						         *
 * 											 *	
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */	

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include "ode_system.h"

/*
 * The function init_y assigns the initial value conditions to the dependent variable 
 * for the "time" t = 0    
 */

int init_y( double t, double y[], void* params )
{
	if ( t == 0 )
	{
		y[0] = 0.5;
		y[1] = 1.0;
		y[2] = -1.0/3.0;
	}
	else
	{
		fprintf( stderr, "Error: initial conditions are only known for t = 0\n" );
		exit( 1 );
	}

	return 0;
}

int ode_system( double t, 
		const double y[], 
		double dydt[],
		void* params )
{
   dydt[0] = 3./5. - 2. *  y[0] - (9./5.) * y[2];
   dydt[1] = 2./5. + 2. *  y[0] - y[1] + (9./5.) * y[2];
   dydt[2] = -2./5.* ( 2 * y[2] + 1 );          

   return GSL_SUCCESS;   
}

char* ode_filename( void )
{
        int name_len = 128;
        char *name = malloc( name_len * sizeof(char) );
        snprintf( name, name_len * sizeof(char), "toy" );
        return name;
}

int write_output( FILE *fp, double t, double t_scale, double y[], void *params )
{
	fprintf( fp, "%.16f %32.16e %32.16e %32.16e \n", t * t_scale, y[0], y[1], y[2] );
	return 0;
}

int ode_term_cond( double t, double y[], void* params )
{
	if ( t > 50.0 )
		return 1;
	else
		return 0;
}

int check_args( void )
{
	return 0;
}
