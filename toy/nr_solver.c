/* 
 *  nr_solver.c 									 
 *  Program numerically evolves a system of ODE's using Numerical Recipes 	 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <nr.h>
#include <nrutil.h>
#include <getopt.h>
#include "ode_system.h"
#include "parse_args.h"

#define N 3	/* number of variables */            
	
/* variables parsed from the command line */
double t_min  = -1;
double t_max  = -1;
double t_step = -1;
double t_scale = 1.0;
double eps    = -1;
double hguess = -1;
char solver_type[SOLVER_TYPE_MAX];
char ode_system_name[ODE_SYSTEM_NAME_MAX];
double e_init = -1;
double f_init = -1;	/* set initial freq. to = 13.3333 */
double m1 = -1;		
double m2 = -1;

int (*init_y)(double, double[], void* );
int (*ode_system)(double, const double[], double dydt[], void* );
int (*ode_term_cond)(double, double y[], void* );
int (*check_args)(void);
char* (*ode_filename)(void);
int (*write_output)( FILE*, double, double, double[], void* );

/* external variable defined in odeint.c  */
float dxsav;	/* smallest interval for which results are stored 		*/
float *xp;  	/* vector of intermediate independent variable values 		*/
float **yp;	/* matrix of intermediate values of the dependedent variables 	*/ 
int kmax;	/* maximum number of steps to be stored				*/
int kount;      /* actual number of steps required to achieve hdid = h		*/  
int nrhs;

void derivs( float nr_t, float nr_y[], float nr_dydt[] )
{

	double t = (double) nr_t;
	double y[3];
	double dydt[3];
	double params[2];
	
	params[0] = m1;
	params[1] = m2;

	y[0] = (double) nr_y[1];
	y[1] = (double) nr_y[2];
	y[2] = (double) nr_y[3];

	ode_system( t, y, dydt, (void*) params );

	nr_dydt[1] = (float) dydt[0];
	nr_dydt[2] = (float) dydt[1];
	nr_dydt[3] = (float) dydt[2];

	nrhs++;
	
	return;
} 

void init_vals( float ystart[], double ode_params[] )
{
	
	double y[3];
	double t = 0.;
	init_y( t, y, (void*) ode_params );

	ystart[1] = (float) y[0];
	ystart[2] = (float) y[1];
	ystart[3] = (float) y[2];

	return;
}

int write_output_nr( FILE *fp, float t, float t_scale, float nr_y[] )
{
	double y[3];	
	double params[2];

	params[0] = m1;
	params[1] = m2;

	y[0] = (double) nr_y[1];
	y[1] = (double) nr_y[2];
	y[2] = (double) nr_y[3];

	return write_output( fp, (double) t, (double) t_scale, y, (void*) params );
}


void ode_termination( float nr_t, float nr_y[] )
{
	double t = (double) nr_t;
	double y[3];

	y[0] = (double) nr_y[1];
	y[1] = (double) nr_y[2];
	y[2] = (double) nr_y[3];

	if ( ode_term_cond ( t, y, NULL ) == 1 )
	{
		fprintf( stderr, " Termination condition reached  \n" );
		exit( 1 );
	}
	
	return;		
}


int main( int argc, char* argv[] )
{
	int nbad; 	/* number of acceptable step tries (hdid == h)   */
	int  nok;	/* number of unacceptable step tries (hdid != h) */ 	
	int i, j; 	/* counters					 */
	int i_max;
	int n;		
	float *ystart;		/* initial values */  
	float *y;		/* stores rows of yp[1][1...kmax], yp[2][1...kmax], yp[3][1...kmax] */ 	
	float t; 		/* time step */
	float t_next;		/* the next time step */ 
	float t_interval; 	/* the difference t_max - t_min */  
	float hmin;		/* min value of integration step */
	float h1;		/* guess for the initial integration step value */
 	double ode_params[4];

	FILE* fp_out = NULL;
#if 0
	FILE* fp_int = NULL;
#endif
	char* of_name;
	const int fname_len = 256;
	char fname[fname_len];

	/* call parse args to get the params we need */ 
	parse_args( argc, argv );

	ode_params[0] =  m1;
	ode_params[1] =  m2;
	ode_params[2] = f_init;
	ode_params[3] = e_init;

	/* create the file name and file pointer */
	of_name = ode_filename();
	snprintf( fname, fname_len * sizeof(char), "nr_solver_%s_%s_%s_%.1e.txt", ode_system_name, of_name, solver_type, eps);
	free( of_name );
	fp_out = fopen( fname, "w" ); 
#if 0
	snprintf( fname, fname_len * sizeof(char), "nr_solver_%s_%s_intermediate_steps_%.1e.txt", ode_system_name, solver_type, eps);
	fp_int = fopen( fname, "w" );
#endif
	
	/* memory allocatation */
	kmax  = 1000000;
	ystart = vector ( 1, N );               
	y      = vector (1, kmax );                   
	xp     = vector ( 1, kmax );                 
	yp     = matrix ( 1, N, 1, kmax );

	nrhs = 0;
	
	/* call init_vals to get the values of initial y values */ 
	t = t_min;
	init_vals( ystart, ode_params );	

	for ( n = 1; n <= N; ++n ) 
	{
  		y[n] = ystart[n]; 
	}	

	hmin = 0.0;
	t_interval = t_max - t_min;
	i_max = (int) ceil( t_interval / t_step );

	/* scale the steps (if requested) */
        t_min /= t_scale;
        t_max /= t_scale;
        t_step /= t_scale;
        if ( hguess > 0 ) h1 = (float) hguess;
        else h1 = t_step;
        fprintf( stderr, "h1 = %32.16e\n", h1 );

	/* print the initial values to file */
#if 0
	fprintf( fp_int, "%.16f %32.16e %32.16e %32.16e \n", t * t_scale, y[1], y[2], y[3] );
#endif
        write_output_nr( fp_out, t, t_scale, y );

	/* check for nans and infs */
	if ( isnan( y[1] ) || isnan( y[2] ) || isnan( y[3] ))
	{
		fprintf( stderr, "y[1] is %f\t, y[2] is %f\t, y[2] is %f\n",
				y[1], y[2], y[3] );
		
	}
	else if ( isinf( y[1] ) || isinf( y[2] ) || isinf( y[3] ))
	{
		fprintf( stderr, "y[1] is %f\t, y[2] is %f\t, y[3] is %f\n",
				y[1], y[2], y[3] );
	}

	/* loop until termination condition is reached */ 
	
	for ( i = 1; i <= i_max; ++i )
	{	
		t        = t_step * (double) i;
		t_next   = t_step * (double) (i + 1);
		dxsav    = (float) ( t_interval / kmax );  
		
		/* call the solver to get the next time step */
		
		if ( ! strcmp( "rkqs", solver_type ) )
			odeint ( ystart, N, (float)t, (float)t_next, (float)eps, 
				h1, hmin, &nok, &nbad, derivs, rkqs );

		else if ( ! strcmp( "bsstep", solver_type ) )
			odeint ( ystart, N, (float)t, (float)t_next, (float)eps, 
				h1, hmin, &nok, &nbad, derivs, bsstep );

		else
		{
			fprintf( stderr, "Error: unknown solver type: %s\n", solver_type );
			exit( 1 );
		}
		
#if 0
		for ( j = 1; j <= kount; ++j )
		{
			/* output to file the stored intermediate values */	
			fprintf( fp_int, "%.16f %32.16e %32.16e %32.16e \n", xp[j] * t_scale, yp[1][j], yp[2][j], yp[3][j] );
		}
#endif

		for ( n = 1; n <= N; n++ )
		{	
			yp[n][i] = ystart[n]; 
			y[n]     = yp[n][i];
		}
#if 0
			fprintf( stderr, "%d sucessful steps\n", nok );
			fprintf( stderr, "%d unsucessful steps\n", nbad );
			fprintf( stderr, "%d function evaluations\n", nrhs );
			fprintf( stderr, "%d stored intermediate values\n", kount );
#endif
			/* call write_output_nr to print the data at uniform time intervals */
                        write_output_nr( fp_out, t, t_scale, y );

			/* call ode_termination to end the integration when it's complete */
			ode_termination( t, y );	
	}

	free_matrix( yp, 1, N, 1, kmax );
	free_vector( xp, 1, kmax);
	free_vector( y, 1, N );
	free_vector( ystart, 1, N );

#if 0
	fclose(fp_int);
#endif
	fclose(fp_out);

	return 0; 
}

