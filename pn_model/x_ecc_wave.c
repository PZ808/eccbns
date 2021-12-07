/* $Id: x_ecc_wave.c,v 1.17 2009/05/20 16:01:41 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "x_ecc_wave.h"

const int num_puts = 1;
const double c_si_pow_4 = C_SI*C_SI*C_SI*C_SI;

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
  double *h_cross_amp_vec;
  double *h_plus_amp_vec;
  double *x_coordinate_vec;
  double *y_coordinate_vec; 

  /* input parameters */
  int con_pn_order;
  int rad_pn_order;
  double mass1, mass2;
  double e_init;
  double f_gw_init;
  double mean_anom_init;
  double ode_eps;
  double sampling_rate;

  /* computed variables */
  double f_gw_isco;
  double omega_init;
  double x_init;
  double x_final, p_final;
  double total_mass;
  double reduced_mass, sym_mass_ratio;

  double amp_plus_cos, amp_plus_sin;
  double amp_cross_cos, amp_cross_sin;

  /* time stepping variables */
  double t, dt, t_min, t_max, t_interval, t_next;
  double step_size;

  /* parameters for the ODE system */
  struct ode_parameters ecc_params;

  /* distance normalization */
  double R = 1.0e-30;
  /* waveform factors */
  double h_factor;    /* units of (time/length)^2 */

  /* iota is inclination angle of orbital plane of *
   * binary with respect to the plane of the ski   */
  const double iota = M_PI/4.0;

  /* parse command line arguments */
  if ( argc != 10 ) 
  {
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, "usage: %s cPN rPN m1 m2 f_init e_init "
        "l_init ode-eps sampling-rate \n", argv[0] );
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

  total_mass = mass1 + mass2;
  reduced_mass = mass1*mass2 / total_mass;
  sym_mass_ratio = reduced_mass / total_mass;

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
      "_%d_%d_%.1f_%.1f_%.2f_%2.1e_%.1f.txt",
      con_pn_order, rad_pn_order, total_mass, f_gw_init,
      e_init, ode_eps, sampling_rate );
  
  strncpy( fname, "x_dyn", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[0] = fopen( fname, "w" );

  strncpy( fname, "x_traj", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[1] = fopen( fname, "w" );

  strncpy( fname, "x_gwave", fname_length );
  strncat( fname, fname_base, fname_length );
  fp[2] = fopen( fname, "w" );

  /* initialize the gsl ode solver */
  double absolute_step_error = 0.0;
  double relative_step_error = ode_eps;
  double y_step_scaling = 1.0;
  double dydt_step_scaling = 1.0;

  const gsl_odeiv_step_type *solver_type = gsl_odeiv_step_rkf45;
  gsl_odeiv_step *solver_step = gsl_odeiv_step_alloc( solver_type, 4 );
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

  /* grab memory for the variables needed for h(t) */
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
  h_plus_amp_vec  = (double*)calloc(i_max, sizeof(double));
  h_cross_amp_vec = (double*)calloc(i_max, sizeof(double));
  
  /* termination gw frequency at isco (geometrized) */
  f_gw_isco = 1.0 / (6.0 * sqrt(6.0) * M_PI * total_mass );
  /* terminate when f = f_isco */
  x_final = pow( M_PI * total_mass * f_gw_isco, 2./3. );
  p_final = 6.0;

  /* initial conditions on dynamical variables and time */
  t_vec[0]   = t    = t_min;
  x_vec[0]   = y[0] = x_init;
  e_vec[0]   = y[1] = e_init;
  l_vec[0]   = y[2] = mean_anom_init;  /* initial mean anomaly */
  phi_vec[0] = y[3] = 0.0;             /* inital phase */
	p_vec[0] = (1.0 - e_vec[0]*e_vec[0]) / pow( M_PI * total_mass * M_SUN_S * f_gw_init, 2./3. );
		
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
		/* p = (1-e*e) / (M*n)^(2/3), where M*n = x^(3/2) in the 0pN RR case */
		p_vec[i] = ( 1.0 - e_vec[i]*e_vec[i] ) / ( x_vec[i] ); 

    /* compute the current value of the eccentric anomaly */
    u_vec[i] = mikkola_finder( e_vec[i], l_vec[i] ); 

    /* compute the values of r using Eq. (5) */
    r_vec[i] = separation( con_pn_order, u_vec[i], sym_mass_ratio, 
        x_vec[i], e_vec[i] );
    
    /* get values for the trajectory */
    x_coordinate_vec[i] = x_cartesian( r_vec[i], phi_vec[i] );
    y_coordinate_vec[i] = y_cartesian( r_vec[i], phi_vec[i] );

    /* check the termination condtion */
    if ( (con_pn_order < 1) ? (p_vec[i] <= 6.0) : (x_vec[i] >= x_final) )
    {
      max_i_reached = i;
			fprintf( stderr, "Terminating at isco: t=%2.2f, x = %2.2f, p=%2.2f\n", 
					t * (total_mass*M_SUN_S), x_vec[i], p_vec[i] );
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
    fprintf( fp[0], 
        "%d %24.16e %24.16e %24.16e %24.16e %24.16e "
        "%24.16e %24.16e %24.16e %24.16e %24.16e\n", 
        (int)i, t_vec[i] * (total_mass*M_SUN_S), x_vec[i], e_vec[i], l_vec[i], 
        phi_vec[i], phi_dot_vec[i], u_vec[i], r_vec[i], r_dot_vec[i], p_vec[i] );

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

    h_factor = ( total_mass * sym_mass_ratio * G_NEWT ) / ( c_si_pow_4 * R );

    /* 
     * The leading order (quadrupolar) post-Newtonian GW polarizations
     * Damour, Gopakumar, and Iyer (PRD 70 064028).
     */

    h_plus_vec[i] = - h_factor * ( 
        (cos(iota)*cos(iota) + 1.0) * ( ( (G_NEWT * total_mass) / r_vec[i] + 
          r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] - 
          r_dot_vec[i] * r_dot_vec[i] ) * cos(2.0 * phi_vec[i]) + 
        2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] * sin(2.0 * phi_vec[i]) ) + 
      ( (G_NEWT * total_mass) / r_vec[i] - r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] 
          - r_dot_vec[i] * r_dot_vec[i] ) * sin(iota)*sin(iota) );

    h_cross_vec[i] = -( 2.0 * h_factor * cos(iota) ) * ( 
        ( (G_NEWT * total_mass) / r_vec[i] + 
         r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] - 
           r_dot_vec[i] * r_dot_vec[i] ) * sin(2.0 * phi_vec[i]) - 
             2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] * cos(2.0 * phi_vec[i])  );

    /* amplitude of the cos(2*phi) term in h_plus */
    amp_plus_cos = -h_factor * ( 1.0 + cos(iota)*cos(iota) ) * ( (G_NEWT * total_mass) / r_vec[i] + 
        r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] - r_dot_vec[i] * r_dot_vec[i] ); 

    /* amplitude of the sin(2*phi) term in h_plus */
    amp_plus_sin = -h_factor * ( 1.0 + cos(iota)*cos(iota) ) * 2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i];

    /* amplitude of the cos(2*phi) term in h_cross */
    amp_cross_cos = -2.0 * h_factor * cos(iota) * (-2.0 * r_vec[i] * r_dot_vec[i] * phi_dot_vec[i] ); 

    /* amplitude of the sin(2*phi) term in h_cross */
    amp_cross_sin = -2.0 * h_factor * cos(iota) * ( (G_NEWT*total_mass)/r_vec[i] + 
        r_vec[i]*r_vec[i]*phi_dot_vec[i]*phi_dot_vec[i] - r_dot_vec[i] * r_dot_vec[i] );
    
    /* amplitude of h_plus */
    h_plus_amp_vec[i] = amp_plus_cos + amp_plus_sin - h_factor * ( 
        (G_NEWT * total_mass) / r_vec[i] - 
       r_vec[i] * r_vec[i] * phi_dot_vec[i] * phi_dot_vec[i] -
        r_dot_vec[i] * r_dot_vec[i] ) * sin(iota)*sin(iota);

    /* amplitude of h_crosss */
    h_cross_amp_vec[i] = amp_cross_cos + amp_cross_sin;
          
    fprintf( fp[2], "%d %24.16e %24.16e %24.16e %24.16e %24.16e\n",
        (int)i, t_vec[i], h_plus_vec[i], h_cross_vec[i], h_plus_amp_vec[i],
        h_cross_amp_vec[i] );
  }

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

  return 0;
}
