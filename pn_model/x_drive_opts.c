/* $Id: x_drive_opts.c,v 1.1 2009/06/02 18:13:24 pzimmerman Exp $ */

/* Computes the gravitational wave signal from 
 * an eccentric compact binary inspiral using 
 * the x-based formalism of Hinder et al. The
 * eccentric waveforms are computed up to and 
 * including 3PN in the conservative dynamics
 * and 2PN in the radiation reaction. 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "x_drive.h"

int parse_args( int argc, char *argv[] );

/* Euler angles (radian  measure)              
 * Converion is that of K.S Thorne (300 Years) 
 */

/* variables from the command line */
int con_pn_order      = -1;
int rad_pn_order      = -1;
double mass1          = -1;
double mass2          = -1;
double e_init         = -1;
double f_gw_init      = -1;
double mean_anom_init = -1;
double ode_eps        = -1;
double sampling_rate  = -1;

const double euler_iota  = 0.0;      /* source's first polar angle       */
const double euler_beta  = 0.0;      /* source's second polar angle      */
const double euler_theta = 0.0;      /* detector's polar angle           */
const double euler_phi   = 0.0;      /* detector's azimuthal             */
const double euler_psi   = 0.0;      /* polarization plane's polar angle */

const double c_si_pow_4 = C_SI*C_SI*C_SI*C_SI;
const int num_puts = 1.0;

int main( int argc, char *argv[] )
{
  unsigned long int i, j, i_max, max_i_reached;
  int badnumber = 0;
  int status = 0;

  /* dynamical quantities */
  double y[4];       /* x, e, l, phi */
  double y_dot[4];   /* x_dot, e_dot, l_dot, phi_dot */

  /* vectors to store the data */
  double *t_vec;
  double *x_vec;
  double *e_vec;
  double *l_vec;
  double *u_vec;
  double *phi_vec;
  double *phi_dot_vec;
  double *r_vec;
  double *r_dot_vec;
	double *p_vec;
  double *h_plus_vec;
  double *h_cross_vec;
  double *h_vec;
  double *x_coordinate_vec;
  double *y_coordinate_vec; 

#if 0
  /* input parameters */
  int con_pn_order;
  int rad_pn_order;
  double mass1, mass2;
  double e_init;
  double f_gw_init;
  double mean_anom_init;
  double ode_eps;         
  double sampling_rate;
#endif

  /* computed variables */
  double f_gw_isco;
  double omega_init;
  double x_init;
  double x_final;
  double x_final_2pn;
  double p_final;
  double num_cycles;
  double total_mass;
  double reduced_mass, sym_mass_ratio;
  double f_plus, f_cross;  /* beam-pattern factors */

  /* time stepping variables */
  double t, dt, t_min, t_max, t_interval, t_next;
  /* step size for the ODE stepping routine */
  double step_size;

  /* parameters for the ODE system */
  struct ode_parameters ecc_params;

  /* distance normalization */
  double R = 1.0e-30;
  /* overall factor in waveform polarizations */
  double h_factor;    /* units are (time/length)^2 */

  /* parse command line arguments */
  parse_args( argc, argv );  
#if 0
  if ( argc != 10 ) 
  {
    fprintf( stderr, "ERROR: Incorrect number of arguments.\n" );
    fprintf( stderr, "Usage: %s\n "
                     "conservative pn order: either 0, 1, 2, or 3\n "
                     "radiation pn order: either 0, 2, 3, or 4\n " 
                     "mass 1 (Msun)\n " 
                     "mass 2 (Msun)\n " 
                     "initial GW frequency (Hz)\n "
                     "initial eccentricity\n "
                     "initial mean-anomaly\t \n "
                     "relative error tolerance of GSL ODE solver\n "
                     "sampling-rate (Hz)\n", 
                     argv[0] );
    return 1;
  }

  con_pn_order   = atoi( argv[1] );
  rad_pn_order   = atoi( argv[2] );
  mass1          = atof( argv[3] );
  mass2          = atof( argv[4] );
  f_gw_init      = atof( argv[5] );
  e_init         = atof( argv[6] );
  mean_anom_init = atof( argv[7] );
  ode_eps        = atof( argv[8] );
  sampling_rate  = atof( argv[9] );
#endif

  /* check the command line arguments */
  if ( e_init < 0.0 || e_init >= 1.0 )
  {
    fprintf( stderr, "Error: invalid arrgument to initial eccentricity\n "
        "eccentricities must be in range [0.0, 1.0)\n" ); 
    exit( 1 );
  }
  if ( (e_init < 0.1) && (ode_eps < 1.0e-16) )
  {
    fprintf( stderr, "Error: invalid argument combination for e_init and ode_eps\n "
        "Either increase e_init above 0.1 with tolerace fixed,\n "
        "or fix e_init and make tolerance looser than 1e-16\n" );
    exit ( 1 );
  }
    
  /* Starting at 1.5 PN radiation reaction we 
   * loosen the tolerance on the ODE solver due to
   * the cost of computing Bessel functions in series. 
   */ 

  if (e_init)
  {
    if (rad_pn_order > 2)
      (e_init <= 0.1) ? (ode_eps *= 1.0e8) : (ode_eps *= 1.0e4);
  }                                                            

  total_mass = mass1 + mass2;
  reduced_mass = mass1*mass2 / total_mass;
  sym_mass_ratio = reduced_mass / total_mass;

  /* numerical factor used in gravitational wave equations */
  h_factor = ( total_mass * sym_mass_ratio * G_NEWT ) /
    ( c_si_pow_4 * 80.0 * R );

  /* Beam-pattern factors for LIGO 
   * Reference: K.S. Thorne "300 Years"
   */
  f_plus = 0.5 * (1.0 + pow( cos(euler_theta), 2.)) * 
    cos(2.* euler_phi) * cos(2.* euler_psi) -
    cos(euler_theta) * sin(2.* euler_phi) * sin(2.* euler_psi);
  f_cross  = 0.5 * (1.0 + pow(cos(euler_theta), 2.)) * 
    cos(2.0*euler_phi) * sin(2.0*euler_psi) +
    cos(euler_theta) * sin(2.0*euler_phi) * cos(2.0*euler_psi);

	/* initial orbital angular frequency in units of solar mass */
	omega_init =  M_PI * f_gw_init * M_SUN_S;  	/* circular orbits f_gw = 2 f_orb */
	x_init = pow( total_mass * omega_init, 2./3. );

  /* store the mass and pn params in the param structure */
  ecc_params.eta = sym_mass_ratio;
  ecc_params.conservative_pn_order = con_pn_order;
  ecc_params.radiation_pn_order = rad_pn_order;

  /* output file pointer and name */
  FILE *fp[3] = {NULL, NULL, NULL};
  const int fname_length = 256;
  char fname_base[fname_length];
  char fname[fname_length];

  snprintf( fname_base, fname_length * sizeof(char),
      "_%d_%d_%.1f_%.1f_%.2f_%2.1e.txt",
      con_pn_order, rad_pn_order, mass1, mass2, e_init, ode_eps );

  strncpy( fname, "xdyn", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[0] = fopen( fname, "w" );

  strncpy( fname, "xtraj", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[1] = fopen( fname, "w" );

  strncpy( fname, "xwave", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[2] = fopen( fname, "w" ); 

  /* GSL ode solver error tolerances */
  double absolute_step_error = 0.0;
  double relative_step_error = ode_eps;
  /* GSL ode solver step scaling */
  double y_step_scaling = 1.0;
  double dydt_step_scaling = 1.0;

  /* GSL Runge-Kutta Fehlberg 4-5 ode stepping method */
  const gsl_odeiv_step_type *solver_type = gsl_odeiv_step_rkf45;
  /* memory for solver */
  gsl_odeiv_step *solver_step = gsl_odeiv_step_alloc( solver_type, 4 );
  /* GSL standard solver control */
  gsl_odeiv_control *solver_control = gsl_odeiv_control_standard_new(
        absolute_step_error, relative_step_error,
        y_step_scaling, dydt_step_scaling );
  gsl_odeiv_evolve *solver_evolve = gsl_odeiv_evolve_alloc( 4 );
  gsl_odeiv_system solver_system = {
    eccentric_x_model_odes, NULL, 4, &ecc_params };
  
  /* set the span of the (maximum) integration time */
  t_min = 0.0;
  t_max = 256.0;

  /* sampling interval [sec] */
  dt = 1.0 / sampling_rate;

  /* compute maximum integration time and maximum time index */
  t_interval = t_max - t_min;
  i_max = (unsigned long int) ceil( t_interval / dt );

  /* scale the start, end and step times by t_sun */
  dt /= total_mass*M_SUN_S;
  t_min /= total_mass*M_SUN_S;
  t_max /= total_mass*M_SUN_S;

  /* grab memory for variables */
  t_vec       = (double*)calloc(i_max, sizeof(double));
  x_vec       = (double*)calloc(i_max, sizeof(double));
  e_vec       = (double*)calloc(i_max, sizeof(double));
  phi_vec     = (double*)calloc(i_max, sizeof(double));
  l_vec       = (double*)calloc(i_max, sizeof(double));
  p_vec       = (double*)calloc(i_max, sizeof(double));
 
  u_vec       = (double*)calloc(i_max, sizeof(double));
  phi_dot_vec = (double*)calloc(i_max, sizeof(double));
  r_vec       = (double*)calloc(i_max, sizeof(double));
  r_dot_vec   = (double*)calloc(i_max, sizeof(double));

  x_coordinate_vec = (double*)calloc(i_max, sizeof(double));
  y_coordinate_vec = (double*)calloc(i_max, sizeof(double));
  
  h_plus_vec  = (double*)calloc(i_max, sizeof(double));
  h_cross_vec = (double*)calloc(i_max, sizeof(double));
  h_vec       = (double*)calloc(i_max, sizeof(double));
  
  /* termination condtions */
  /* gw frequency (geometrized) */
  f_gw_isco = 1.0 / (6.0 * sqrt(6.0) * M_PI * total_mass );
    

  /* If 1pN (RR) or less we terminate when x is at x(f=f_gw_isco) = 1/6 */
  x_final = pow( M_PI * total_mass * f_gw_isco, 2./3. );
  /* If 2pN (RR) we must use the 2PN isco */
  x_final_2pn = 3.0*sym_mass_ratio * (1.0 - sqrt(1.0 - (14.0*sym_mass_ratio)/9.0));

  /* If the dynamics are "Newtonian" we terminate 
   * when p is at p(f=f_gw_isco) = 6 */
  p_final = 6.0;

  /* initial conditions on dynamical variables and time */
  t_vec[0]   = t    = t_min;
  x_vec[0]   = y[0] = x_init;
  e_vec[0]   = y[1] = e_init;
  l_vec[0]   = y[2] = mean_anom_init;  /* initial mean anomaly */
  phi_vec[0] = y[3] = 0.0;             /* inital phase */
	p_vec[0] = (1.0 - e_vec[0]*e_vec[0]) / 
    pow( M_PI * total_mass * M_SUN_S * f_gw_init, 2./3. );
		
  /* initial eccentric anomaly */
  u_vec[0] = pn_kepler_equation( con_pn_order, sym_mass_ratio, x_vec[0], e_vec[0], l_vec[0] );

  /* compute the intital value of r using Eq. (5) */
  r_vec[0] = separation( con_pn_order, u_vec[0], sym_mass_ratio, 
      x_vec[0], e_vec[0] );

  /* compute the initial values for the trajectory */
  x_coordinate_vec[0] = x_cartesian(r_vec[0], phi_vec[0]);
  y_coordinate_vec[0] = y_cartesian(r_vec[0], phi_vec[0]);
  
  /* evolve the dynamical variables forward until we reach  */
  /* termination condition: either t > t_max or x > x_final */
  for ( i = 1; i < i_max; ++i )
  {
    /* set the current time, next time and step size for the ode solver */
    t = dt * (double) (i-1);
    t_next = dt * (double) i;
    step_size = dt;

    /* integrate to the next time step */
    while ( t < t_next )
    {
      status = gsl_odeiv_evolve_apply(
        solver_evolve, solver_control,
        solver_step, &solver_system,
        &t, t_next, &step_size, y );

      /* check the the integrator worked */
      if ( status != GSL_SUCCESS )
      {
        fprintf( stderr, "failure in GSL integrator\n" );
        return 1;
      }
    }
    /* check for nan or inf in dynamical variables */
    for ( j = 0; j < 4; ++j )
    {
      if ( isnan( y[j] ) || isinf( y[j] ) )
      {
        badnumber = 1;
        max_i_reached = i;
      }
    }

    /* stop integrating if we got a nan or an inf */
    if ( badnumber ) 
		{
			fprintf( stderr, "Terminating due to nan of inf at t = %2.2f, x = %2.2f, p = %2.2f, x_final = %2.2f\n", 
					t * (total_mass*M_SUN_S), x_vec[i], p_vec[i], x_final );
			break;
		}
    
    /* compute derivatives at the current value of t */
    eccentric_x_model_odes( t, y, y_dot, (void*) &ecc_params );

    /* store the computed time, dynamical variables and their derivatives */
    t_vec[i] = t;
    x_vec[i] = y[0];
    e_vec[i] = y[1];
    l_vec[i] = y[2];
    phi_vec[i] = y[3];
    phi_dot_vec[i] = y_dot[3];

		/* semi-latus rectum: correct only for 0pN R.R. 0pN Con. */
		/* p = (1-e*e) / (M*n)^(2/3), where M*n = x^(3/2)  */
		p_vec[i] = ( 1.0 - e_vec[i]*e_vec[i] ) / ( x_vec[i] ); 

    /* compute the current value of the eccentric anomaly */
    u_vec[i] = pn_kepler_equation( con_pn_order, sym_mass_ratio,
        x_vec[i], e_vec[i], l_vec[i] ); 

    /* compute the values of r using Eq. (5) */
    r_vec[i] = separation( con_pn_order, u_vec[i], sym_mass_ratio, 
        x_vec[i], e_vec[i] );
    
    /* get values for the trajectory */
    x_coordinate_vec[i] = x_cartesian( r_vec[i], phi_vec[i ]);
    y_coordinate_vec[i] = y_cartesian( r_vec[i], phi_vec[i] );

    /* check the termination condtion for low rad-PN case */
    if ( (rad_pn_order <= 2) && ((con_pn_order < 1) ? (p_vec[i] <= 6.0) : (x_vec[i] >= x_final)) )
    {
      max_i_reached = i;
#if 0
      fprintf( stderr, "Terminating at 1PN-isco: t = %2.2f, x = %.2f, p = %2.2f\n", 
          t * (total_mass*M_SUN_S), x_vec[i], p_vec[i] );
#endif
      num_cycles = phi_vec[i]/M_PI;
      fprintf( stderr, "Terminating at schwarzschild-isco: t = %2.2f, x = %.2f, GW_cycles = %g\n", 
          t * (total_mass*M_SUN_S), x_vec[i], num_cycles );
      break;
    }

     /* Inner most stable orbits for rad PN orders greater than 1PN 
      * http://arxiv.org/abs/gr-qc/0209089v2 
      * 1PN ISCO: x = 1/6
      * 2PN ISCO: x = 3/(14*eta) * ( 1-sqrt(1 - (14*eta)/9 )
      * taylor expanding the 2PN result gives 
      * x = (1/6)*( 1 + (7/18)*eta + Order(eta*eta) )
      */

    if ( (rad_pn_order > 2) && (x_vec[i] >= x_final) )
    {
      max_i_reached = i;
      num_cycles = phi_vec[i]/M_PI;
      fprintf( stderr, "Terminating at schwarzschild-isco: t = %2.2f, x = %.2f, GW_cycles = %g\n", 
          t * (total_mass*M_SUN_S), x_vec[i], num_cycles );
      break;
    }

  } /* end integration for loop */

  /* compute the numerical derivative of r(t) */
  for ( i = 2; i < max_i_reached-2; ++i )
  {
    r_dot_vec[i] = 
      ( -r_vec[i+2] + 8.0*r_vec[i+1] - 8.0*r_vec[i-1] + r_vec[i-2] ) /
      ( 12.0 * dt );
  }  

  /* store the intermediate results in an output file */
  for ( i = 0; i < max_i_reached; ++i )
  {
    /* print the dynamical variables to file */
    fprintf( fp[0], 
        "%d %24.16e %24.16e %24.16e %24.16e %24.16e "
        "%24.16e %24.16e %24.16e %24.16e %24.16e\n", 
        (int)i, t_vec[i] * (total_mass*M_SUN_S), x_vec[i], e_vec[i], l_vec[i], 
        phi_vec[i], phi_dot_vec[i], u_vec[i], r_vec[i], r_dot_vec[i], p_vec[i] );

    /* print the trajectory data to file */
    if ( (i % num_puts) == 0 )
    {
      fprintf( fp[1], "%24.16e %24.16e %24.16e\n", 
        t_vec[i] * (total_mass*M_SUN_S), x_coordinate_vec[i], y_coordinate_vec[i] );
      fflush( fp[1] );
    }
  }

  /* compute the gravitational wave signal */
  for ( i = 2; i < max_i_reached-2; ++i )
  {
    /* convert to SI units */
    t_vec[i] *= (total_mass*M_SUN_S);
    r_vec[i] *= (total_mass*M_SUN_M);
    r_dot_vec[i] *= (M_SUN_M/ M_SUN_S);
    phi_dot_vec[i] /= (total_mass*M_SUN_S);


    /* 
     * The leading order (quadrupolar) post-Newtonian GW polarizations
     * Reference: Damour, Gopakumar, and Iyer (PRD 70 064028).
     */

    h_plus_vec[i] = - h_factor * ( 
        (cos(euler_iota)*cos(euler_iota) + 1.0) * ( ( (G_NEWT * total_mass) / r_vec[i] + 
          r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] - 
          r_dot_vec[i] * r_dot_vec[i] ) * cos(2.0 * (phi_vec[i]-euler_beta)) + 
        2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] * sin(2.0 * phi_vec[i]) ) + 
      ( (G_NEWT * total_mass) / r_vec[i] - r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] 
          - r_dot_vec[i] * r_dot_vec[i] ) * sin(euler_iota)*sin(euler_iota) );

    h_cross_vec[i] = -( 2.0 * h_factor * cos(euler_iota) ) * ( 
        ( (G_NEWT * total_mass) / r_vec[i] + 
         r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] - 
           r_dot_vec[i] * r_dot_vec[i] ) * sin(2.0*(phi_vec[i]-euler_beta)) - 
             2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] * cos(2.0*(phi_vec[i]-euler_beta))  );

    /* Compute gravitational wave signal */
    h_vec[i] = f_plus*h_plus_vec[i] + f_cross*h_cross_vec[i];
                         
    /* Print waveform data to file */
    fprintf( fp[2], "%d %24.16e %24.16e %24.16e %24.16e\n",
        (int)i, t_vec[i], h_plus_vec[i], h_cross_vec[i], h_vec[i] );
  }

  /* Close the files and free memory */
  for( j = 0; j < 3; ++j )
  {
    fclose( fp[j] );
  }

  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);

  free( t_vec );
  free( x_vec );
  free( e_vec );
  free( phi_vec );
  free( l_vec );
 
  free( u_vec );
  free( phi_dot_vec );
  free( r_vec );
  free( r_dot_vec );

  free( x_coordinate_vec );
  free( y_coordinate_vec );
  
  free( h_plus_vec );
  free( h_cross_vec );
  free( h_vec );

  return 0;
}
