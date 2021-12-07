#include <stdio.h>
#include <math.h>
#include <qm.h>
#include <string.h>
#include <qmtypes.h>
#include <except.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>       
#include "t4.h"

int t4_ode_system( double t, const double y[], 
    double dydt[], void *params )
{
  const double euler_gamma = 0.57721566490153286;

  struct t4_params *tp = (struct t4_params *) params;

  double M   = tp->mass1 + tp->mass2;
  double mu  = tp->mass1 * tp->mass2 / M;
  double eta = mu/M;
  double x   = y[0];
  double phi = y[1];

  double x_pow_5 = x*x*x*x*x;

  /*  \dot{x} in geometrized units  */
  dydt[0] = ( 64.0*eta*x_pow_5 / (5.0*M) ) * 
    ( 1.0 - (743.0/336.0 + 11.0*eta/4.0) * x +
             (4.0*QM_PI) * x*sqrt(x)
  (34103.0/18144.0 + 13661.0*eta/2016.8 + 
   59.0*eta*eta/18.0) *x*x );

  /* \dot{\phi} in geometrized units */
  dydt[1] = x*sqrt(x)/M;

  return GSL_SUCCESS;
}                             

int t4_waveform( real_vec_t** h_plus, real_vec_t** h_cross, 
                  double m1,
                  double m2,
                  double f_min, 
                  int N, 
                  double dt 
                  )    
{

  const int numvars = 2;
  int badnumber = 0;
  int i, j, max_i_reached;
  int errcode = QM_NO_ERROR; 

  /* time stepping variables */
  double t, t_next;

  /* step size returned by the ODE stepping routine */
  double step_size;

  double m_tot, mu, eta;
  double amp_factor;      

  /* dynamical evolution variables and their derivatives */
  double y[numvars], dydt[numvars];
  double f_gw_init, f_gw_isco;
  double omega_init;
  double x, x_init, x_final, x_final_2pn;
  double phi;

  /* set the mass vaiables */
  m_tot = m1+m2;
  mu = m1*m2 / m_tot;
  eta = mu / m_tot;

  /* scale the step size by the total mass in units of the sun in seconds */
  const double geometrized_m_tot = m_tot * QM_MTSUN_SI;
  /* cannonical distance for waveform */
  const double dist_factor = -4.0 * mu * QM_MRSUN_SI * DYN_RANGE_FAC 
    / (1.0e6 * QM_PC_SI); 
  /* time step in units of seconds / geometrized total mass */
  const double step = dt / geometrized_m_tot;
  
  /* parameters for the ODE system */
  struct t4_params ode_params;

  ode_params.mass1 = m1;                                        
  ode_params.mass2 = m1;                                        
  
/* GSL ode solver error tolerances */
  double absolute_step_error = 0.0;
  double relative_step_error = 1.0e-6;
  /* GSL ode solver step scaling */
  double y_step_scaling = 1.0;
  double dydt_step_scaling = 1.0;    

  /* GSL Runge-Kutta Fehlberg 4-5 ode stepping method */
  const gsl_odeiv_step_type *solver_type = gsl_odeiv_step_rkf45;
  gsl_odeiv_step *solver_step = gsl_odeiv_step_alloc( solver_type, numvars ); 
  gsl_odeiv_control *solver_control = gsl_odeiv_control_standard_new( 
      absolute_step_error, relative_step_error, y_step_scaling, dydt_step_scaling );    
  gsl_odeiv_evolve *solver_evolve = gsl_odeiv_evolve_alloc( numvars );
  gsl_odeiv_system solver_system = { t4_ode_system, NULL, numvars, &ode_params };

  /* If there is no memory preallocated for the output data then create a  */
  /* vector of length N. Otherwise check that the existing memory is of    */
  /* length N. This allows the function to be called from python, where    */
  /* the memory needs to be allocated, or from a C function which has      */
  /* pre-allocated the memory. If the memory has been preallocated, set it */
  /* it to zero to avoid it conting garbage.                               */
  if ( ! (*h_plus) )
    *h_plus = new_real_vec_t( NULL, N );
  else if ( (*h_plus)->n == N )
    memset( (*h_plus)->d, 0, N * sizeof(float) );
  else
    errcode = QM_VALUE_ERROR;

  if ( ! (*h_cross) )
    *h_cross = new_real_vec_t( NULL, N );
  else if ( (*h_cross)->n == N )
    memset( (*h_cross)->d, 0, N * sizeof(float) );
  else
    errcode = QM_VALUE_ERROR;

  /* check that all the input arrays are the correct length */
  if ( errcode == QM_VALUE_ERROR )
  {
    throw_exception( QM_VALUE_ERROR, "output vectors are incorrect length" );
    return QM_VALUE_ERROR;
  }

  /* store the sample interval in the output vectors */
  (*h_plus)->dx = (*h_cross)->dx = dt;

  /* initial conditions */
  /* initial orbital angular frequency in units of solar mass */
  omega_init =  QM_PI * f_gw_init * QM_MTSUN_SI;  	
  x_init = pow( m_tot * omega_init, 2./3. );            
  /* initial conditions on dynamical variables */
  y[0] = x_init;  
  y[1] = phi = 0.0;  

  f_gw_isco = 1.0 / (6.0 * sqrt(6.0) * QM_PI * m_tot );
  /* final value of x is at x(f_isco) = 1/6 */
  x_final = pow( QM_PI * m_tot * f_gw_isco, 2./3. );
  /* 2PN isco if we need it (diffenent from schwarschild) */
  x_final_2pn = 3.0 * eta * 
    (1.0 - sqrt(1.0 - 14.0 * eta / 9.0)); 
  
  t = 0;
  i = 0;

 /* We evolve the dynamical variables forward until we reach the
  * termination condition: x > x_final                     
  */    

  while( 1 )
  {
    /* if we have run out of memory for the waveform, break out of the loop */
    if ( i >= N  )
    {
      throw_exception( QM_MEMORY_ERROR, "output too short for waveform" );
      errcode = QM_MEMORY_ERROR;
      break;
    }
    /* compute the gravitational waveform from the dynamical */
    /* variables and store it in the output structures       */
     amp_factor = x * dist_factor;

    (*h_plus)->d[i]  = (float) ( amp_factor * cos(2.0*phi) );
    (*h_cross)->d[i] = (float) ( amp_factor * sin(2.0*phi) );
    
    /* advance the time */
    t = i * step;
    t_next = ++i * step;
    step_size = step;

    /* call the solver to get the next timestep */
    while ( t < t_next && errcode == GSL_SUCCESS )
      errcode = gsl_odeiv_evolve_apply(
          solver_evolve, solver_control, solver_step, &solver_system,
          &t, t_next, &step_size, y );
    
    /* check for that the solver exited successfully */
    if ( errcode != GSL_SUCCESS )
    {
      throw_exception( QM_RUNTIME_ERROR, "integrator failure" );
      errcode = QM_RUNTIME_ERROR;
      break;
    }

    /* copy the output variables necessary to compute the gw */
    x   = y[0];
    phi = y[1];
    
    /* check for nan or inf in dynamical variables */
    for ( j = 0; j < 3; ++j )
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
      throw_exception( QM_RUNTIME_ERROR, "NaN in dynamical variable" );
      errcode = QM_RUNTIME_ERROR;
      break;
    }
    /*  terminate at ISCO */
    if ( x >= x_final )
      break;

  } /* end of integration loop */

  gsl_odeiv_evolve_free (solver_evolve);
  gsl_odeiv_control_free (solver_control);
  gsl_odeiv_step_free (solver_step);        

  return QM_NO_ERROR;
}
