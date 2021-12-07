/* $Id: eccentric_orbit_newt.c,v 1.8 2009/02/14 22:28:43 pzimmerman Exp $ */

/*
 * 	eccentric_orbit_newt.c                               
 * 	Numerically evolves Mercury
 */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "kepler.h"

#define T_SUN 4.92549095e-6
#define DAY 86400.0

const double AU = 1.49597870691e11;

/* mass of sun in solar units */
double mass_sun = 1.0;
/* mass of mercury in solar units */
double mass_merc = 1.7634e-7;
double eccentricity = 0.2056;
/* semi-major axis in SI units */
double semi_major_axis = 5.791e10;     
/* period of mercury in SI units */
double period = 7.6005438e6;

double mikkola_finder( double e, double l);

double relative_separation( double ecc_anomaly )
{
  return semi_major_axis * ( 1.0 - eccentricity * cos( ecc_anomaly ) );
}

double relative_dphi_dt( double mean_motion,
  	double ecc_anomaly )
{
  return mean_motion * sqrt( 1.0 - eccentricity*eccentricity ) /
  	pow( 1.0 - eccentricity * cos( ecc_anomaly ), 2.0 );
}

/* Cartesian coordinates (x, y) are 
 * taken with the origin at the center 
 * of the ellipse.
 * When we introduce pN corrections we  
 * should  measure the coordinates 
 * relative to the center of mass (the focus).
 */

double x_cartesian( double ecc_anomaly )
{
	return semi_major_axis * cos(ecc_anomaly);    
}

double y_cartesian( double ecc_anomaly )
{
	double semi_minor_axis = semi_major_axis * sqrt(1.0 - eccentricity*eccentricity);
	return semi_minor_axis * sin(ecc_anomaly);
}

/* define the ode system */
int keplerian_evolution ( double t, const double y[],
		double dydt[], void *params )
{ 
	double *kepler_parameters = (double*)params;
  double mean_motion  = kepler_parameters[0];
  double ecc_anomaly  = kepler_parameters[1];

  dydt[0] = relative_dphi_dt(mean_motion, ecc_anomaly);

	return GSL_SUCCESS;
}

int main( int argc, char *argv[] )
{
	char *root_finder;
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int out_interval;
	int status;

	double y[1];
	double kepler_parameters[2];
	double x_coordinate, y_coordinate;
	double r,	x_au, y_au;
	double P_1, P_2, dP_1, dP_2;
	double total_mass;
	double mean_anomaly, mean_motion, ecc_anomaly;
	double mean_anomaly_init;
	double tolerance;
	double ode_eps;
	double t, dt, t_min, t_max, t_interval, t_next;
	double t_days, out_days;
	double step, step_size;
	double sampling_rate, sampling_interval;

	total_mass = mass_sun + mass_merc;

	if ( argc != 9 )
  {
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, 
     "usage: %s t-min t-max initial-mean-anomaly tolerance ode-eps sampling-rate output-interval root-finder \n", argv[0] );
    return 1;
  }

	t_min 			  		 = atof( argv[1] );
	t_max 			  		 = atof( argv[2] );
	mean_anomaly_init  = atof( argv[3] );
	tolerance 	  		 = atof( argv[4] );
	ode_eps            = atof( argv[5] );
	sampling_rate      = atof( argv[6] );
	out_interval       = atoi( argv[7] ); 
	root_finder        = argv[8];

	t_days = t_max/DAY;
	out_days = out_interval/DAY;
	mean_motion = (2.0 * M_PI * T_SUN) / period;
	kepler_parameters[0] = mean_motion;

	/* output file parameters */
	FILE *fp = NULL;
	const int fname_length = 256;
	char fname[fname_length];

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
			"zeroPN_mercury_evolution_%s_%.2f_%2.1e_%2.1e_%.1f_%g.txt", 
			root_finder, t_days, tolerance, ode_eps, sampling_rate, (float)out_interval);
	fp = fopen (fname, "w");

	double absolute_step_error = 0.0;
	double relative_step_error = ode_eps;
	double y_step_scaling = 1.0;
	double dydt_step_scaling = 1.0;

	/* sampling interval [sec] */
	sampling_interval = 1.0 / sampling_rate;

	/* scale by mass of the sun in seconds */
	dt = sampling_interval/T_SUN;
  step = dt;

	/* GSL differential equation solving mechanism */
	const gsl_odeiv_step_type *solver_type
		= gsl_odeiv_step_rkf45;
 	gsl_odeiv_step *solver_step
		= gsl_odeiv_step_alloc( solver_type, 1 );
	gsl_odeiv_control *solver_control
		= gsl_odeiv_control_standard_new( 
				absolute_step_error, relative_step_error,
				y_step_scaling, dydt_step_scaling );
	gsl_odeiv_evolve *solver_evolve
		= gsl_odeiv_evolve_alloc( 1 );
	gsl_odeiv_system solver_system = { 
		keplerian_evolution, NULL, 1, 
		(void*)kepler_parameters };

	t = t_min;
	t_min /= T_SUN;
	y[0] = 0.;
	t_max /= T_SUN;
  t_interval = t_max - t_min;
  i_max = (unsigned long int) ceil( t_interval / dt ); 

	/* step the time in units of M_SUN */
	for (i = 1; i <= i_max; i++)
	{
		t = step * (double)i;
		t_next = step * (double)(i+1);
		step_size = step;

		mean_anomaly = mean_anomaly_init + mean_motion * t;

  	if ( !strcmp("guide", root_finder) ) 
  		ecc_anomaly = kepler_guide(eccentricity, mean_anomaly, tolerance);
  	else if ( !strcmp("mikkola", root_finder) )
   		ecc_anomaly = mikkola_finder(eccentricity, mean_anomaly);
  	else
  	{
    	fprintf( stderr, "Error: unknown solver type: %s\n", root_finder );
     	exit( 1 );
  	}
 
  	r = relative_separation(ecc_anomaly);
		x_coordinate = x_cartesian(ecc_anomaly);
		y_coordinate = y_cartesian(ecc_anomaly);
		/* compute period from mean_anomaly */
		P_1 = (2.0 * M_PI * t) / mean_anomaly;
		/* compute period from dphi/dt */
		P_2 = ( 2.0 * M_PI * sqrt(1.0 - eccentricity*eccentricity) ) /
		 ( (1.0 - eccentricity*cos(ecc_anomaly))*(1.0 - eccentricity*cos(ecc_anomaly)) * 
		  relative_dphi_dt(mean_motion, ecc_anomaly) );

		/* check for nan of inf */
    if ( isnan(y[0]) || isinf(y[0]) )
		{
			fprintf( stderr, "y[0] is %f\n", y[0] );
			break;
		}

		kepler_parameters[1] = ecc_anomaly;

 		status = gsl_odeiv_evolve_apply(
    	solver_evolve, solver_control,
    	solver_step, &solver_system,
    	&t, t_next, &step_size, y );

    if ( status != GSL_SUCCESS )
    {
    	fprintf( stderr, "failure in GSL integrator\n" );
    	break;
    }

    if ( i % out_interval == 0 ) 
    {
    	x_au = x_coordinate/AU; 
    	y_au = y_coordinate/AU; 
    	r /= AU;
    	t *= T_SUN;
    	P_1 *= T_SUN;
    	P_2 *= T_SUN;
    	dP_1 = fabs(P_1-period) / DAY;
    	dP_2 = fabs(P_2-period) / DAY;

			fprintf( fp, "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %12.12e %12.12e\n", 
				t/DAY, y[0], r, x_au, y_au, ecc_anomaly, mean_anomaly, dP_1 , dP_2);
		}
	}	

	fclose( fp );
  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);
	return 0;
}
