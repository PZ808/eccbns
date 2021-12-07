/* $Id: x_con.c,v 1.2 2009/04/02 13:43:28 pzimmerman Exp $ */

/*
 *	evolves post-Newtonian conservative dynamics only
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "conservative_ode_system.h"
#include "conservative_dynamics_pure.h"
#include "x_con.h"

#define DAY 86400.0
/* mass of the sun in kg */
#define M_SUN 1.9884e30
/* mass of the sun in seconds */
#define T_SUN 4.92549095e-6 
/* mass of the sun in meters */
#define L_SUN 1.47662504e3
/* Newton's gravitational constant */
#define G_NEWT 6.6743e-11

const double AU = 1.49597870691e11;


double eccentricity = 0.2;

/* semi-major axis in SI units */
double semi_major_axis = 1.0e6;     

int conservative_evolution_equations ( double t, const double y[],
		double dydt[], void *params )
{ 
	double *kepler_parameters = (double*)params;
  double x = kepler_parameters[0];
  double ecc_anomaly = kepler_parameters[1];
  double eta = kepler_parameters[2];
                      
  dydt[0] = dphi_dt(ecc_anomaly, eta, eccentricity, x);

	return GSL_SUCCESS;
}

int main( int argc, char *argv[] )
{
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int out_interval;
	int status;

	double y[1];
	double kepler_parameters[3];
	double mass_1, mass_2;
	double mean_motion;
	double period;
	double precession_angle;
	double reduced_mass, sym_mass_ratio;
	double omega_init, mean_anom_init;
	double x;
	double total_mass;
	double mean_anomaly, ecc_anomaly;
	double ode_eps;
	double r;
	double t, dt, t_min, t_max, t_interval, t_next, t_days;
	double step, step_size, sampling_rate, sampling_interval;

	if ( argc != 9 ) {
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, 
    		"usage: %s m1 m2 t-min t-max initial-mean-anomaly ode-eps sampling-rate output-interval\n", argv[0] );
    return 1;
  }

	mass_1			   = atof( argv[1] );
	mass_2			   = atof( argv[2] );
	t_min 			   = atof( argv[3] );
	t_max 			   = atof( argv[4] );
	mean_anom_init = atof( argv[5] );
	ode_eps        = atof( argv[6] );
	sampling_rate  = atof( argv[7] );
	out_interval   = atoi( argv[8] ); 
  
	total_mass = mass_1 + mass_2;
	reduced_mass = mass_1*mass_2 / total_mass;
	sym_mass_ratio = reduced_mass / total_mass;

	kepler_parameters[2] = sym_mass_ratio;

	/* orbital period in seconds from Kepler's Law*/
	period = 2.0*M_PI * sqrt( 
	 (semi_major_axis*semi_major_axis*semi_major_axis) /
			(G_NEWT*total_mass*M_SUN) );

  /* x and e are constant in the conservative dynamics */
  omega_init = ( (2.0 * M_PI * T_SUN) / period );
  x = pow(total_mass * omega_init, 2./3.);
  kepler_parameters[0] = x;

	/* mass scaled mean_motion */
  mean_motion = mean_motion_func(eccentricity, sym_mass_ratio, x);

  /* precession angle per orbit Eqn. 4.15 D&D */
  precession_angle = ( (6.0 * M_PI * L_SUN * total_mass) /
  		(semi_major_axis * (1.0 - eccentricity*eccentricity)) ); 

	/* output file parameters */
	FILE *fp = NULL;
	const int fname_length = 256;
	char fname[fname_length];
	t_days = t_max/DAY;

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
		"cpn_x_model_evolution_%.1f_%.1f_%.1f_%g_%.2fsec_%2.1e_%.1f_%g.txt", 
		mass_1, mass_2, eccentricity, semi_major_axis, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp = fopen (fname, "w");

	double absolute_step_error = 0.0;
	double relative_step_error = ode_eps;
	double y_step_scaling = 1.0;
	double dydt_step_scaling = 1.0;

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
		conservative_evolution_equations, NULL, 1, 
		(void*)kepler_parameters };

	y[0] = 0.0;
	/* sampling rate [Hz] */
	sampling_interval = 1.0 / sampling_rate;
	/* scale the steps by t_sun */
	dt = sampling_interval/T_SUN;
  step = dt;
	t = t_min;
	t_min /= T_SUN;
	t_max /= T_SUN;
  t_interval = t_max - t_min;
  i_max = (unsigned long int) ceil( t_interval / dt );
	
	for (i = 1; i <= i_max; i++)
	{
		t = step * (double)i;
		t_next = step * (double)(i+1);
		step_size = step;

		mean_anomaly = mean_anom_init + mean_motion * t;
		ecc_anomaly = mikkola_finder( eccentricity, 
				mean_anomaly );
		r = separation(ecc_anomaly, sym_mass_ratio, eccentricity, x); 

		kepler_parameters[1] = ecc_anomaly;

 		status = gsl_odeiv_evolve_apply(
      solver_evolve, solver_control,
     	solver_step, &solver_system,
     	&t, t_next, &step_size, y );

		/* check for nan of inf */
    if ( isnan(y[0]) || isinf(y[0]) )
		{
			fprintf( stderr, "y[0] is %f\n", y[0] );
			break;
		}
    if ( status != GSL_SUCCESS )
    {
     	fprintf( stderr, "failure in GSL integrator\n" );
     	break;
    }
    if ( i % out_interval == 0 ) 
    {
    	//r *= L_SUN; 
    	t *= T_SUN;
    	double phase = y[0];
    	double num_cycles = (phase/(2.0*M_PI + precession_angle));
    	double x_pha = r * cos(phase);
    	double y_pha = r * sin(phase); 

			fprintf( fp, 
				"%12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n", 
				t, phase, r, x_pha, y_pha, ecc_anomaly, 
				num_cycles, num_cycles*precession_angle, num_cycles*(M_PI/30.));
			fflush( fp );
		}
	}	

	fclose( fp );
  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);
	return 0;
}
