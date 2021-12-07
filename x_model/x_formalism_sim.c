/* $Id: x_formalism_sim.c,v 1.7 2009/04/09 23:04:38 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "x_formalism_sim.h"
#include "conservative_dynamics.h"
#include "reactive_dynamics.h"

/* mass of the sun in kg */
#define M_SUN 1.9884e30
/* mass of the sun in seconds */
#define T_SUN 4.92549095e-6 
/* mass of the sun in meters */
#define L_SUN 1.47662504e3
/* Newton's gravitational constant */
#define G_NEWT 6.6743e-11
/* speed of light */
#define C_SI 299792458.0


/* function: dphi_dt
 * computes  angular frequency 
 */
double dphi_dt(
		int conservative_pn_order,
		double ecc_anomaly, 
		double eta,
		double x, 
		double eccentricity )
{
	double phi_dot_0pn = rel_dphi_dt_0pn(
			eccentricity,
			ecc_anomaly,
			eta );
	double phi_dot_1pn = rel_dphi_dt_1pn(
			eccentricity,
			ecc_anomaly,
			eta );
	double phi_dot_2pn = rel_dphi_dt_2pn(
			eccentricity,
			ecc_anomaly,
			eta );
	/* mass scaled dphi/dt */
  if (conservative_pn_order == 0)
		return( phi_dot_0pn * x*sqrt(x) );
  else if (conservative_pn_order == 1)
		return ( ( phi_dot_0pn + phi_dot_1pn * x) * x*sqrt(x) );   
  else if (conservative_pn_order == 2)
		return ( ( phi_dot_0pn + phi_dot_1pn * x + phi_dot_2pn * x*x ) * x*sqrt(x) );   
  else 
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit(1);
	}   
}

/* function: separation 
 * computes orbital separation 
 */

double separation( 
		int conservative_pn_order,
		double ecc_anomaly, 
		double eta,
		double x,
		double eccentricity)
{
	double r_0pn = rel_sep_0pn( 
			eccentricity, 
			ecc_anomaly );
	double r_1pn = rel_sep_1pn( 
			eccentricity, 
			ecc_anomaly, 
			eta ); 
	double r_2pn = rel_sep_2pn( 
			eccentricity, 
			ecc_anomaly, 
			eta ); 
#if 0
	double r_3pn = rel_sep_3pn( 
			eccentricity, 
			ecc_anomaly, 
			eta ); 
#endif

  if (conservative_pn_order == 0)
		return( r_0pn );
  else if (conservative_pn_order == 1)
		return ( r_0pn * (1./x) + r_1pn );
  else if (conservative_pn_order == 2)
		return ( r_0pn * (1./x) + r_1pn + r_2pn * x );
  else 
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit(1);
	}   
}

double mean_motion( 
		int conservative_pn_order,
		double eta,
		double x,
		double eccentricity ) 
{	
	/* mass scaled mean_motion */
  double n_0pn = pow(x, 3./2.);
	double n_1pn = mean_motion_1pn( eccentricity );
	double n_2pn = mean_motion_2pn( eccentricity, eta );

  if (conservative_pn_order == 0)
		return( n_0pn );
  else if (conservative_pn_order == 1)
		return ( (1.0 + n_1pn*x) * sqrt(x)*x );
  else if (conservative_pn_order == 2)
		return ( (1.0 + n_1pn*x + n_2pn*x*x) * sqrt(x)*x );
  else 
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit(1);
	}   
}

double x_cartesian( double r, double phi )
{
	return (r * cos(phi));    
}

double y_cartesian( double r, double phi )
{
	return (r * sin(phi));
}

int eccentric_x_model_odes( double t, const double y[],
		double dydt[], void *params )
{ 
	struct ode_parameters *op = (struct ode_parameters *)params;
	double eta = op->eta;
  double ecc_anomaly = op->ecc_anomaly;
  int conservative_pn_order = op->conservative_pn_order;
  int radiation_pn_order = op->radiation_pn_order;

  /* mass scaled radiative ODEs 
   * dxdt = M*x_dot (A25)
   * dedt = M*e_dot (A30)
   */
  double dxdt;
  double dedt;
  /* mass scaled radiative ODEs 
   * dldt = M*l_dot (8)
   * dphidt = M*phi_dot (A10)
   */ 
	double dldt;
	double dphidt;

  double x = y[0];
  double e = y[1];

  //double l   = y[2];
  //double phi = y[3];

  if ( radiation_pn_order == 0 )
	{
		dxdt = x_dot_0pn(e,eta) * x*x*x*x*x;
		dedt = e_dot_0pn(e,eta) * x*x*x*x;
	}
	else if ( radiation_pn_order == 1 )
	{
		dxdt = ( x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x ) * x*x*x*x*x;
		dedt = ( e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x ) * x*x*x*x*x;
	}
	else if ( radiation_pn_order == 2 )
	{
		dxdt = ( x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x 
				+ x_dot_2pn(e, eta) * x*x ) * x*x*x*x*x;
		dedt = ( e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x 
				+ e_dot_2pn(e, eta) * x*x ) * x*x*x*x*x;
	}
	else 
	{
		fprintf( stderr, "Error in PN order: %d\n", radiation_pn_order );
		exit(1);
	}

	dldt = mean_motion( conservative_pn_order, eta, x, e);
	dphidt = dphi_dt(conservative_pn_order, ecc_anomaly, eta, x, e);

	dydt[0] = dxdt;
	dydt[1] = dedt;
	dydt[2] = dldt;
  dydt[3] = dphidt;

	return GSL_SUCCESS;
}

int main( int argc, char *argv[] )
{
	int j;
	int badnumber;
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int max_count;
	int status;

	/* 
	 * dynamical quantities 
	 */

	double y[4];
	double y_dot[4];
	/* vectors to store the data */
	double *t_vec;
	double *phi_vec;
	double *phi_dot_vec;
	double *r_vec;
	double *r_dot_vec;
	double *ell_vec;
	/* variable for kepler's equation */
	double mean_anomaly, ecc_anomaly;
	double sep, phase;
	/* phase space postion coordinates */ 
	double x_coordinate, y_coordinate; 

	/* 
	 * input parameters 
	 */

	int con_pn_order;
	int rad_pn_order;
	double mass1, mass2;
	double ode_eps;
	double f_gw_init;  
  double mean_anom_init;
	double e_init;
  double sampling_rate;
	unsigned long int out_interval;

	/* computed variables */
	double total_mass;
	double reduced_mass, sym_mass_ratio;
	double omega_init;
	double x_init; 
	double f_gw_isco; 
	double x_final;
	/* time stepping variables */
	double t, dt, t_min, t_max, t_interval, t_next;
	double step, step_size, sampling_interval;

  /* waveforms */
  double h_plus, h_cross;

  /* distance normalization */
  double R = 1.0e-30;
  /* inclination angle */
  const double iota = M_PI/4.0;

	/* parameters for the ODE system */
	struct ode_parameters ecc_params;

  /* parse the command line arguments */
	if ( argc != 11 ) 
	{
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, 
    		"usage: %s cPN rPN m1 m2 f_init e_init l_init ode-eps samp-rate out-interval\n", argv[0] );
    return 1;
  }

	con_pn_order   = atoi( argv[1] );
	rad_pn_order   = atoi( argv[2] );
	mass1 			   = atof( argv[3] );
	mass2	  		   = atof( argv[4] );
	f_gw_init 		 = atof( argv[5] );
	e_init 	 			 = atof( argv[6] );
	mean_anom_init = atof( argv[7] );
	ode_eps        = atof( argv[8] );
	sampling_rate  = atof( argv[9] );
	out_interval   = atoi( argv[10] ); 

	total_mass = mass1 + mass2;
	reduced_mass = mass1*mass2 / total_mass;
	sym_mass_ratio = reduced_mass / total_mass;

	 /* initial omega in terms of M-sun */ 
  omega_init = ( 0.5 * M_PI * f_gw_init * T_SUN );
  /* initial value of x (geometrized) */
  x_init = pow(total_mass * omega_init, 2./3.);
	/* duration of simulation in seconds */
	t_max = 256.0;
	/* gw frequency at isco (geometrized) */
	f_gw_isco = (1.0 / (6.0 * sqrt(6.0) * total_mass ) );
  /* terminate when f = f_isco */
  x_final = pow( 0.5 * M_PI * total_mass * f_gw_isco, 2./3. );
	fprintf( stderr, "x_final = %e\n", x_final );

	ecc_params.eta = sym_mass_ratio;
	ecc_params.conservative_pn_order = con_pn_order;
	ecc_params.radiation_pn_order = rad_pn_order;

	/* output file parameters */
	FILE *fp = NULL;
	FILE *fp_traj  = NULL;
	FILE *fp_dat_a = NULL;
	FILE *fp_dat_b = NULL;
	const int fname_length = 256;
	char fname[fname_length];

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
			"x_dyn_%d_%d_%.1f_%.1f_%.1f_%.2f_%2.1e_%.1f_%g.txt", 
			con_pn_order, rad_pn_order, total_mass, f_gw_init, 
			e_init, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp = fopen (fname, "w");

	snprintf(fname, fname_length * sizeof(char),
			"traj_%d_%d_%.1f_%.1f_%.1f_%.2f_%2.1e_%.1f_%g.txt", 
			con_pn_order, rad_pn_order, total_mass, f_gw_init, 
			e_init, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp_traj = fopen (fname, "w");

	snprintf(fname, fname_length * sizeof(char),
			"storedvals_%d_%d_%.1f_%.1f_%.1f_%.2f_%2.1e_%.1f_%g.dat", 
			con_pn_order, rad_pn_order, total_mass, f_gw_init, 
			e_init, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp_dat_a = fopen (fname, "w");

	snprintf(fname, fname_length * sizeof(char),
			"x_wave_%d_%d_%.1f_%.1f_%.1f_%.2f_%2.1e_%.1f_%g.dat", 
			con_pn_order, rad_pn_order, total_mass, f_gw_init, 
			e_init, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp_dat_b = fopen (fname, "w");

	double absolute_step_error = 0.0;
	double relative_step_error = ode_eps;
	double y_step_scaling = 1.0;
	double dydt_step_scaling = 1.0;

	/* GSL differential equation solving mechanism */
	const gsl_odeiv_step_type *solver_type
		= gsl_odeiv_step_rkf45;
 	gsl_odeiv_step *solver_step
		= gsl_odeiv_step_alloc( solver_type, 4 );
	gsl_odeiv_control *solver_control
		= gsl_odeiv_control_standard_new( 
				absolute_step_error, relative_step_error,
				y_step_scaling, dydt_step_scaling );
	gsl_odeiv_evolve *solver_evolve
		= gsl_odeiv_evolve_alloc( 4 );
	gsl_odeiv_system solver_system = { 
		eccentric_x_model_odes, NULL, 4, 
		&ecc_params };

  /* initial conditions on y */
	y[0] = x_init;
	y[1] = e_init;
	y[2] = mean_anom_init;	/* initial mean anomaly */
	y[3] = 0.0;	  					/* inital phase */

  /* initial eccentric anomaly */
	ecc_anomaly = mikkola_finder( e_init, mean_anom_init );

	/* sampling interval [sec] */
	sampling_interval = 1.0 / sampling_rate;
	/* scale the steps by t_sun */
	dt = sampling_interval/T_SUN;
  step = dt;
	t = t_min;
	t_min /= T_SUN;
	t_max /= T_SUN;
  t_interval = t_max - t_min;
  i_max = (unsigned long int) ceil( t_interval / dt );

	//fprintf( fp_dat_b, "%% dx = %16.8e", dt );

  /* grab memory for the variables needed for h(t) */
  t_vec 	  	= (double*)calloc(i_max, sizeof(double));
  phi_vec 	  = (double*)calloc(i_max, sizeof(double));
  phi_dot_vec = (double*)calloc(i_max, sizeof(double));
  r_vec 			= (double*)calloc(i_max, sizeof(double));
 	r_dot_vec 	= (double*)calloc(i_max, sizeof(double));
  ell_vec 	 	= (double*)calloc(i_max, sizeof(double));

	/* set the initial conditions of the vectors */
	t_vec[0] = t;
  phi_vec[0] = y[3];
	/* mass scaled initial orbital frequency */
  phi_dot_vec[0] = 0.0;
	/* mass scaled initial relative separation */
	r_vec[0] = separation( con_pn_order, ecc_anomaly, sym_mass_ratio, x_init, e_init ); 
	/* mass scaled initial relative velocity */
	r_dot_vec[0] = 0.0;
	ell_vec[0] = mean_anom_init;

  /* inital x_dot */
	y_dot[0] = 0.0;
  /* inital e_dot */
	y_dot[1] = 0.0;
  /* intial l_dot = n_init = pi * f_gw_init */
	y_dot[2] = M_PI * f_gw_init;
	/* initial phi_dot */
	y_dot[3] = 0.0; 


	for (i = 1; i < i_max; i++)
	{
		t = step * (double)i;
		t_next = step * (double)(i+1);
		step_size = step;

		/* store the mass scaled time */
		t_vec[i] = t;

 		status = gsl_odeiv_evolve_apply(
      	solver_evolve, solver_control,
     		solver_step, &solver_system,
     		&t, t_next, &step_size, y );

    if ( status != GSL_SUCCESS )
    {
     	fprintf( stderr, "failure in GSL integrator\n" );
     	break;
    }

		mean_anomaly = y[2];
		ell_vec[i] = mean_anomaly;
		/* call the root finder to get ecc_anom */
		ecc_anomaly = mikkola_finder( y[1], 
				mean_anomaly );
  	ecc_params.ecc_anomaly = ecc_anomaly;
		/* compute the separation */
		sep = separation( con_pn_order, ecc_anomaly, 
				sym_mass_ratio, y[0], y[1] ); 
		r_vec[i] = sep;
		phi_vec[i] = phase = y[3];

		/* get values for the trajectory */
		x_coordinate = x_cartesian(sep, phase);
		y_coordinate = y_cartesian(sep, phase);

  	/* termination condtion */
  	if ( y[0] == x_final )
  		break;
		
		/* check for nan of inf */
		for ( j = 0; j < 4; ++j ) 
    	if ( isnan( y[j] ) || isinf( y[j] ) )
			{
				badnumber = 1;
				max_count = i;
			}
    if ( badnumber ) 
      break;

		//fprintf( stderr, "y[%d] is %f\n", j, y[j] );
		/* call eccentric_x_model_odes to compute 
		 * the derivatives dxdt, dedt, dldt, and dphidt */
		eccentric_x_model_odes( t, y, y_dot, (void*)&ecc_params );
    phi_dot_vec[i] = y_dot[3];

		/* write output */
    if ( i % out_interval == 0 ) 
    {
			t *= T_SUN;
			fprintf( fp, "%24.16e %24.16e %24.16e %24.16e %24.16e\n", 
					t, y[0], y[1], y[2], y[3] );

			fprintf( fp_traj, "%24.16e %24.16e\n", x_coordinate, y_coordinate );
		}
	}	

	//i_max = max_count;

	/* compute dr/dt using the 2-point algorithm */
	for ( i = 1; i < i_max-1; ++i )
	{
		r_dot_vec[i] = ( r_vec[i+1] - r_vec[i-1] ) / ( 2.0 * dt );

    fprintf( fp_dat_a, "%d %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
				i, t_vec[i]*T_SUN, r_vec[i], r_dot_vec[i], phi_vec[i], phi_dot_vec[i], ell_vec[i] );
	}

	/* compute h_plus and h_cross */
	for ( i = 1; i < i_max-1; ++i )
	{
		/* convert to SI units */
		r_vec[i] *= L_SUN; 
		r_dot_vec[i] *= (L_SUN/T_SUN);
		phi_dot_vec[i] /= (T_SUN);

		 h_plus	= - ( (total_mass * sym_mass_ratio * G_NEWT) / (C_SI*C_SI*C_SI*C_SI * R)) * 
			( (cos(iota)*cos(iota) + 1.0) *
				( cos( phi_vec[i] ) * ( -r_dot_vec[i]*r_dot_vec[i] + 
				r_vec[i]*r_vec[i] * phi_dot_vec[i]*phi_dot_vec[i] + total_mass/r_vec[i] ) + 
				2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] * sin( 2.* phi_vec[i] ) ) -  
					( -r_dot_vec[i] * r_dot_vec[i] - 
						r_vec[i]*r_vec[i] * phi_dot_vec[i]*phi_dot_vec[i] + 
							total_mass/r_vec[i] ) * sin(iota)*sin(iota) );  

		h_cross	= - ( (2.* G_NEWT * total_mass * sym_mass_ratio) / ( C_SI*C_SI*C_SI*C_SI * R ) ) * cos(iota) *
      ( ( -r_dot_vec[i]*r_dot_vec[i] * phi_dot_vec[i]*phi_dot_vec[i] + 
      		total_mass / r_vec[i] ) * sin( 2.* phi_vec[i] ) -
      	- 2.0 * r_vec[i] * cos(2.* phi_vec[i] ) * r_dot_vec[i] * phi_dot_vec[i] );

    fprintf( fp_dat_b, "%24.16e %24.16e %24.16e\n", t_vec[i]*T_SUN, h_plus, h_cross );
  }

	fclose( fp );
	fclose( fp_traj );
	fclose( fp_dat_a );
	fclose( fp_dat_b );
  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);
	return 0;
}
