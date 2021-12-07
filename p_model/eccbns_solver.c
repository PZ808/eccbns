/* 
 *  eccbns_solver.c
 *  Program numerically evolves a system of ODE's using GSL  	 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include "ode_system.h"
#include "parse_args.h"

/* variables parsed from the command line */
double t_min = -1;
double t_max = -1;
double t_step = -1;
double t_scale = 1.0;	
double eps = -1; 		
double hguess = -1;
char solver_type[SOLVER_TYPE_MAX]; 
char ode_system_name[ODE_SYSTEM_NAME_MAX];
double e_init = -1;
double f_init = -1;	
double m1 = -1;		
double m2 = -1;

int (*init_y)(double, double[], void* );
int (*ode_system)(double, const double[], double dydt[], void* );
int (*ode_term_cond)(double, double y[], void* );
int (*check_args)(void);
char* (*ode_filename)(void);
int (*write_output)( FILE*, double, double, double[], void* );

int main( int argc, char* argv[] )
{
	int status; 
  unsigned long int i, i_max;
	double t, t_next; 
	double t_interval;
	double h;

	/* FILE* fp_test = NULL; */
	FILE* fp_out = NULL;
	char* of_name;
	const int fname_len = 256;
	char fname[fname_len];
	double y[3];
	double ode_parameters[4];

	gsl_odeiv_step_type* T = NULL;

	/* call parse args to get the params we need */ 
	parse_args( argc, argv );

	if ( ! strcmp( "rkf45", solver_type ) )
		T = (gsl_odeiv_step_type*) gsl_odeiv_step_rkf45;
	else if ( ! strcmp( "rkck", solver_type ) )
		T = (gsl_odeiv_step_type*)  gsl_odeiv_step_rkck;
	else
	{
		fprintf( stderr, "Error: unknown solver type: %s\n", solver_type );
		exit( 1 );
	}

	ode_parameters[0] = m1; 
	ode_parameters[1] = m2;
	ode_parameters[2] = f_init;
	ode_parameters[3] = e_init;

	/* create the ODE solver */
	gsl_odeiv_step* solver_step
		= gsl_odeiv_step_alloc(T, 3);
	gsl_odeiv_control* solver_control
		= gsl_odeiv_control_standard_new(0.0, eps, 1.0, 1.0);
	gsl_odeiv_evolve* solver_evolve
		= gsl_odeiv_evolve_alloc (3);
	gsl_odeiv_system solver_system 
		= { ode_system, NULL, 3, (void*) ode_parameters };

	/* create the file name and file pointer */
	of_name = ode_filename();
	snprintf( fname, fname_len*sizeof(char), "pwave_%s_%.1e.txt", of_name, eps);
	free( of_name );
	fp_out = fopen( fname, "w" ); 
#if 0
	of_name = ode_filename();
	snprintf( fname, fname_len*sizeof(char), "p_test_%s_%.1e_%.1f.txt", of_name, eps, 1./t_step );
	free( of_name );
	fp_test = fopen( fname, "w" ); 
#endif

	t = t_min;

	/* call init_y to get the intitial values */ 
	init_y( t, y, (void*) ode_parameters );

	t_interval = t_max - t_min;
	i_max = (unsigned long int) ceil ( t_interval / t_step );

	/* scale the steps (if requested) */
	t_min /= t_scale;
	t_max /= t_scale;
	t_step /= t_scale;

	/* call write output to print the initial values to file */
	write_output( fp_out, t, t_scale, y, (void*) ode_parameters );
	/* fprintf( fp_test, "%24.16e %24.16e %24.16e\n", t*t_scale, y[0], y[2]); */

	/* loop until termination condition is reached */ 
	for( i = 0; i <= i_max ; i++)
	{
		t = t_step * (double) i;
		t_next = t_step * (double) (i + 1);
		if ( hguess > 0 ) h = hguess;
		else h = t_step; 

		/* Check for not a number in dynamical variables */
		if ( isnan( y[0] ) || isnan( y[1] ) || isnan( y[2] ))
		{
			fprintf( stderr, "y[0] is %f\t, y[1] is %f\t, y[2] is %f\n",
					y[0], y[1], y[2] );
			break;
		}
		while ( t < t_next )
		{

			/* call the solver to get the next time step */
			status = gsl_odeiv_evolve_apply ( solver_evolve, solver_control, solver_step, 
					&solver_system,
					&t, t_next,
					&h, y );

			/* check for a succesful solver exit */
			if ( status != GSL_SUCCESS )
			{
				fprintf( stderr, "Error: failure in GSL integrator \n" );
				exit( 1 );
			}     
		}

		/* call write_output to print the data at uniform time intervals */
		write_output( fp_out, t, t_scale, y, (void*) ode_parameters );
    /* print out t, phi, and e to compare with xmodel */
#if 0
		fprintf( fp_test, "%24.16e %24.16e %24.16e\n", t*t_scale, y[0], y[2]);
#endif
		/* call ode_termination to end the integration when it's complete */
		if ( ode_term_cond( t, y, NULL ) )
		{
			fprintf( stderr, "Termination condition reached at t = %f\n", t * t_scale );

			break; 
		}
	}

	fclose( fp_out ); 

	gsl_odeiv_evolve_free( solver_evolve );
	gsl_odeiv_control_free( solver_control );
	gsl_odeiv_step_free( solver_step );

	return 0;
}
