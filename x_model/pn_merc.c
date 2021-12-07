/* $Id: pn_merc.c,v 1.5 2009/03/01 22:54:23 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "conservative_dynamics.h"

#define DAY 86400.0
#define YEAR 31557600.
#define M_SUN 1.98892e30
/* mass of the sun in seconds */
#define T_SUN 4.92549095e-6 
/* mass of the sun in meters */
#define L_SUN 1.47662504e3

const double AU = 1.49597870691e11;

/* mass of sun in solar units */
double mass_sun = 1.0;
/* mass of mercury in solar units */
double mass_merc = 1.7634e-7;
double mu = 1.7634e-7;
double eta = 1.7634e-7;
double eccentricity = 0.2056;
/* semi-major axis in SI units */
double semi_major_axis = 5.791e10;     
/* period of mercury in SI units */
double period = 7.600521e6;

double mikkola_finder( double e, double l);

double dphi_dt(double ecc_anomaly, double x)
{
	double phi_dot_0pn = rel_dphi_dt_0pn( 
			eccentricity, 
			ecc_anomaly, 
			eta );
	double phi_dot_1pn = rel_dphi_dt_1pn( 
			eccentricity, 
			ecc_anomaly, 
			eta );

	/* mass scaled dphi/dt */
	double M_times_phi_dot = phi_dot_0pn * pow(x, 3./2.)
					+ phi_dot_1pn * pow(x, 5./2.);

	return M_times_phi_dot;
}

double separation(double ecc_anomaly, double x)
{
	double r_0pn = rel_sep_0pn(
			eccentricity, 
			ecc_anomaly );
	double r_1pn = rel_sep_1pn(
			eccentricity, 
			ecc_anomaly, 
			eta ); 
	double r_over_M = r_0pn * (1./x) + r_1pn;
	return r_over_M;
}

double x_cartesian( double ecc_anomaly )
{
	return semi_major_axis * cos(ecc_anomaly);    
}

double y_cartesian( double ecc_anomaly )
{
	double semi_minor_axis = semi_major_axis * 
		sqrt(1.0 - eccentricity*eccentricity);
	return semi_minor_axis * sin(ecc_anomaly);
}

int keplerian_evolution ( double t, const double y[],
		double dydt[], void *params )
{ 
	double *kepler_parameters = (double*)params;
  double x = kepler_parameters[0];
  double ecc_anomaly  = kepler_parameters[1];

  dydt[0] = dphi_dt(ecc_anomaly, x);

	return GSL_SUCCESS;
}

int main( int argc, char *argv[] )
{
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int out_interval;
	int status;

	double y[1];
	double kepler_parameters[2];
	double x_coordinate, y_coordinate, x_au, y_au;
	double omega_init, mean_anom_init;
	double x;
	double mean_motion;
	//double P_1, P_2, dP_1, dP_2;
	double total_mass;
	double mean_anomaly, ecc_anomaly;
	double ode_eps;
	double r;
	double n_0pn, n_1pn;
	double t, dt, t_min, t_max, t_interval, t_next, t_days, t_years;
	double step, step_size, sampling_rate, sampling_interval;

	if ( argc != 7 ) {
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, "usage: %s t-min t-max initial-mean-anomaly ode-eps sampling-rate output-interval\n", argv[0] );
    return 1;
  }

	t_min 			   = atof( argv[1] );
	t_max 			   = atof( argv[2] );
	mean_anom_init = atof( argv[3] );
	ode_eps        = atof( argv[4] );
	sampling_rate  = atof( argv[5] );
	out_interval   = atoi( argv[6] ); 
  
	total_mass = mass_sun + mass_merc;
  /* x and e are constant in the conservative dynamics */
  omega_init = ( (2.*M_PI*T_SUN) / period );
  x = pow( total_mass * omega_init, 2./3. );
  kepler_parameters[0] = x;
	/* mean_motion is constant in the conservative dynamics */
  n_0pn = pow(x, 3./2.);
	n_1pn = mean_motion_1pn(eccentricity);
	/* mass scaled mean_motion */
  mean_motion = n_0pn + n_1pn * pow(x, 5./2.);

	/* output file parameters */
	FILE *fp = NULL;
	const int fname_length = 256;
	char fname[fname_length];
	t_days = t_max/DAY;
	t_years = t_max/YEAR;

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
		"pn_mercury_evolution_%.2fyrs_%2.1e_%.1f_%g.txt", 
		t_years, ode_eps, sampling_rate, (float)out_interval);
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
		keplerian_evolution, NULL, 1, 
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
		r = separation(ecc_anomaly, x); 
		x_coordinate = x_cartesian( ecc_anomaly );
		y_coordinate = y_cartesian( ecc_anomaly );

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
    	x_coordinate *= L_SUN;
    	y_coordinate *= L_SUN;
    	x_au = x_coordinate/AU;
    	y_au = y_coordinate/AU;
    	r = (L_SUN/AU) * r; 
    	t = (T_SUN/DAY) * t;
    	double phase = y[0];

			fprintf( fp, 
				"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
				t, phase, r, x_au, y_au, mean_anomaly, ecc_anomaly);
		}
	}	

	fclose( fp );
  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);
	return 0;
}
