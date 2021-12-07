/* $Id: x_formalism.c,v 1.5 2009/04/01 12:14:51 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_roots.h>

#include "x_formalism.h"
#include "conservative_dynamics.h"
//#include "reactive_dynamics.h"
//#include "mean_motion_init.h"

/* mass of the sun in kg */
#define M_SUN 1.9884e30
/* mass of the sun in seconds */
#define T_SUN 4.92549095e-6 
/* mass of the sun in meters */
#define L_SUN 1.47662504e3
/* Newton's gravitational constant */
#define G_NEWT 6.6743e-11

struct ic_params
  {
    double f_init;    
    double M;         
    double eta;     
    double e;       
    int pn_order;
  };

double n_init( double n, void *params );


double x_dot_0pn( double, double );
double x_dot_1pn( double, double );
double x_dot_2pn( double, double );
double e_dot_0pn( double, double );
double e_dot_1pn( double, double );
double e_dot_2pn( double, double );

double dphi_dt(
		int conservative_pn_order,
		double ecc_anomaly, 
		double eta,
		double x, 
		double eccentricity )
{


  double e = eccentricity;
  double u_factor = e*cos(ecc_anomaly) - 1.0;

                                                  
  double phi_dot_0pn = ( sqrt(1.-e*e) / ((u_factor*u_factor)) ); 

  double phi_dot_1pn = ( - (e*(eta - 4.) * u_factor) /
		(sqrt(1. - e*e) * u_factor*u_factor*u_factor) );



	/* mass scaled dphi/dt */
  if (conservative_pn_order == 0)
		return( phi_dot_0pn * x*sqrt(x) );
  else if (conservative_pn_order == 1)
		return ( ( phi_dot_0pn + phi_dot_1pn * x) * x*sqrt(x) );   
 #if 0
  else if (conservative_pn_order == 2)
		return ( ( phi_dot_0pn + phi_dot_1pn * x + phi_dot_2pn * x*x ) * x*sqrt(x) );   
#endif
  else 
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit(1);
	}   
}

/* function: separation 
 * computes orbital separation 
 */
#if 0
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
	double r_3pn = rel_sep_3pn( 
			eccentricity, 
			ecc_anomaly, 
			eta ); 

  /* mass scaled separation */
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
#endif

double mean_motion( 
		int conservative_pn_order,
		double eta,
		double x,
		double eccentricity ) 
{	

	double ecc = eccentricity;
	/* mass scaled mean_motion */
  double n_0pn = pow(x, 3./2.);

	double n_1pn = (3./(ecc-1.0));
  double n_2pn = ( ((26.*eta - 51.)*ecc*ecc + 28. - 18.)
			/ (4.*(ecc*ecc-1.) * (ecc*ecc-1.)) );

  /* mass scaled mean motion*/
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

  double x   = y[0];
  double e   = y[1];

  //double l   = y[2];
  //double phi = y[3];

  if ( radiation_pn_order == 0 )
	{
		dxdt = x_dot_0pn(e,eta)*x*x*x*x*x;
		dedt = e_dot_0pn(e,eta)* x*x*x*x;
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
		fprintf( stderr, "Error in pN order: %d\n", radiation_pn_order );
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
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int out_interval;
	int status;
	int con_pn_order;
	int rad_pn_order;

	double y[4];
	//double kepler_parameters[3];
	double mass1, mass2;
	double total_mass;
	double reduced_mass, sym_mass_ratio;
	double x_coordinate, y_coordinate; 
	double period;
	double f_gw_init;
	double omega_init, mean_anom_init;
	double x_init, e_init;
	double mean_anomaly, ecc_anomaly;
	double ode_eps;
	double r;
	/* time stepping variables */
	double t, dt, t_min, t_max, t_interval, t_next;
	double step, step_size, sampling_rate, sampling_interval;


	if ( argc != 11 ) 
	{
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, 
    		"usage: %s cPN rPN m1 m2 f_init e_init l_init ode-eps samp-rate out-interval\n", argv[0] );
    return 1;
  }

  /* command line arguments */
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

	/* parameters for the ODE system */
	struct ode_parameters ecc_params;

	ecc_params.eta = sym_mass_ratio;
	ecc_params.conservative_pn_order = con_pn_order;
	ecc_params.radiation_pn_order = rad_pn_order;

  /* root solving mechanism 
   * for initial conditions 
   * on the mean motion "n"
   */

  int n_init_status;
  unsigned int iter = 0;
  unsigned int max_iter = 1000;
  double n = 0.0;           	/* starting value (guess) */
  double n_lo = 100.0;        /* estimated lower bracket */
  double n_hi  = 200.0;       /* estimated upper bracket */
  double root_eps_abs = 0.0;
  double root_eps_rel = 1.0e-6;

  const gsl_root_fsolver_type *root_finder_type;
  gsl_root_fsolver *root_solver;

  gsl_function F;
  struct ic_params n_params = { f_gw_init, total_mass, 
  	sym_mass_ratio, e_init, rad_pn_order };
  F.function = &n_init;
  F.params   = &n_params;

  root_finder_type = gsl_root_fsolver_brent;
  root_solver = gsl_root_fsolver_alloc (root_finder_type);
  gsl_root_fsolver_set( root_solver, &F, n_lo, n_hi );

  iter = 0;

  do {
  	n_init_status = gsl_root_fsolver_iterate (root_solver);
  	n = gsl_root_fsolver_root (root_solver);
  	n_lo = gsl_root_fsolver_x_lower (root_solver);
  	n_hi = gsl_root_fsolver_x_upper (root_solver);
  	n_init_status = gsl_root_test_interval ( n_lo, n_hi,
  	    root_eps_abs, root_eps_rel);

  	if ( n_init_status == GSL_SUCCESS )
  	{
  	  fprintf( stderr, "n_initial = %e \n", n*T_SUN );
  	  break;
  	}

  	if ( iter == max_iter || n_init_status != GSL_CONTINUE )
  	{
  	  fprintf( stderr, " Error in GSL root finder \n" );
  	  exit(1);
  	}
  } while (1);
	
	/* orbital period in seconds from pN Kepler's Law */
	period = 2.0*M_PI / n;
  /* initial value of x obtained from n_init */
  omega_init = ( (2.*M_PI*T_SUN) / period );
  x_init = pow(total_mass * omega_init, 2./3.);

	/* output file parameters */
	FILE *fp = NULL;
	const int fname_length = 256;
	char fname[fname_length];

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
			"x_model_waveforms_%.1f_%.1f_%.1f_%.2f_%2.1e_%.1f_%g.txt", 
			total_mass, f_gw_init, e_init, t_max, ode_eps, sampling_rate, (float)out_interval);
	fp = fopen (fname, "w");

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

  /*
   * initial conditions
   */

	t = 0.0;

	y[0] = x_init;
	y[1] = e_init;
	y[2] = n * t;	/* mean anomaly */
	y[3] = 0.0;	  /* phase        */

	/* termination condition 
   * 	???
   */

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

		mean_anomaly = y[2];
		/* call the root finder to get ecc_anom */
		ecc_anomaly = mikkola_finder( y[1], 
				mean_anomaly );
  	
  	ecc_params.ecc_anomaly = ecc_anomaly;

		//r = separation( con_pn_order, ecc_anomaly, 
			//	sym_mass_ratio, y[0], y[1] ); 

    /* put numerical differentiation for dr/dt scheme here */

		x_coordinate = x_cartesian(r, y[3]);
		y_coordinate = y_cartesian(r, y[3]);


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
    	t *= T_SUN; 
			fprintf( fp, 
					"%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
					t, y[0], y[1], y[2], y[3], x_coordinate, y_coordinate, ecc_anomaly);
		}
	}	

	fclose( fp );
  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);
	return 0;
}
