/* $Id: x_ecc_ode.c,v 1.13 2009/05/20 20:57:54 pzimmerman Exp $ */
/* x_ecc_ode.c 
 *  Older version of x_ecc_ode_newt.c 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

#include "x_ecc_wave.h"

int eccentric_x_model_odes(
    double t, 
    const double y[],
    double dydt[], 
    void *params )
{ 
  /* parse the paramater structure */
  struct ode_parameters *op = (struct ode_parameters *)params;
  double eta = op->eta; 
  int conservative_pn_order = op->conservative_pn_order;
  int radiation_pn_order = op->radiation_pn_order;
  
  /* input variables */
  double x = y[0];
  double e = y[1];
  double l = y[2];

  double u = mikkola_finder( e, l );

  dydt[0] = dx_dt( radiation_pn_order, eta, x, e );
  dydt[1] = de_dt( radiation_pn_order, eta, x, e );
  dydt[2] = dl_dt( conservative_pn_order, eta, x, e );

  /* account for the zero eccentricitiy case */
  if ( e )
    dydt[3] = dphi_dt( conservative_pn_order, u, eta, x, e );
  else 
    dydt[3] = phi_dot_0pn(e, eta, u);

  return GSL_SUCCESS;
}

double dx_dt( int radiation_pn_order, double eta, double x, double e )
{
  double x_pow_5 = x*x*x*x*x;
  double xdot = 0;

  if ( radiation_pn_order == 0 ) /* 0 pN term */
  {
    xdot = x_dot_0pn(e, eta) * x_pow_5;
  }
  else if ( radiation_pn_order == 1 ) /* 0.5 pN term */
  {
    xdot = x_dot_0pn(e, eta) * x_pow_5;
  }
  else if ( radiation_pn_order == 2 ) /* 1 pN term */
  {
    xdot = ( x_dot_0pn(e, eta) + x_dot_1pn(e, eta) * x ) * x_pow_5;
  }
  else if ( radiation_pn_order == 3 ) /* 1.5 pN term */
  {
    xdot = ( x_dot_0pn(e,eta) + x_dot_1pn(e, eta) * x + 
        x_dot_1p5pn(e, eta) * x * sqrt(x) ) * x_pow_5;
  }
  else if ( radiation_pn_order == 4 ) /* 2 pN term */
  {
    xdot = ( x_dot_0pn(e,eta) + x_dot_1pn(e, eta) * x + 
        x_dot_1p5pn(e, eta) * x * sqrt(x) +
        x_dot_2pn(e, eta) * x * x ) * x_pow_5;
  }
  else 
  {
    fprintf( stderr, "Error in PN order: %d\n", radiation_pn_order );
    exit(1);
  }

  return xdot;
}

double de_dt( int radiation_pn_order, double eta, double x, double e )
{
  double x_pow_4 = x*x*x*x;
  double edot = 0;

  if ( radiation_pn_order == 0 ) /* 0 pN term */
  {
    edot = e_dot_0pn(e, eta) * x_pow_4;
  }
  else if ( radiation_pn_order == 1 ) /* 0.5 pN term */
  {
    edot = e_dot_0pn(e, eta) * x_pow_4;
  }
  else if ( radiation_pn_order == 2 ) /* 1 pN term */
  {
    edot = ( e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x ) * x_pow_4;
  }
  else if ( radiation_pn_order == 3 ) /* 1.5 pN term */
  {
    edot = ( e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
        e_dot_1p5pn(e, eta) * x * sqrt(x) ) * x_pow_4;
  }
  else if ( radiation_pn_order == 4 ) /* 2 pN term */
  {
    edot = ( e_dot_0pn(e, eta) + e_dot_1pn(e, eta) * x +
        e_dot_1p5pn(e, eta) * x * sqrt(x) + 
        e_dot_2pn(e, eta) * x * x ) * x_pow_4;
  }
  else 
  {
    fprintf( stderr, "Error in PN order: %d\n", radiation_pn_order );
    exit(1);
  }

  return edot;
}

double dl_dt( int conservative_pn_order, double eta, double x, double e ) 
{
  double x_pow_3_2 = sqrt(x) * x;
  double ldot = 0;

  if (conservative_pn_order == 0) /* 0 pN term */
  {
    ldot = x_pow_3_2;
  }
  else if (conservative_pn_order == 1) /* 1 pN term */
  {
    ldot = (1.0 + x * l_dot_1pn(e, eta)) * x_pow_3_2;
  }
  else if (conservative_pn_order == 2) /* 2 pN term */
  {
    ldot = (1.0 + x * l_dot_1pn(e, eta) + 
        x * x * l_dot_2pn(e, eta)) * x_pow_3_2;
  }
  else if (conservative_pn_order == 3) /* 3 pN term */
  {
    ldot = (1.0 + x * l_dot_1pn(e, eta) + 
        x * x * l_dot_2pn(e, eta) +
        x * x * x * l_dot_3pn(e, eta) ) * x_pow_3_2;
  }
  else 
  {
    fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
    exit(1);
  }   

  return ldot;
}

double dphi_dt( int conservative_pn_order, 
    double u, double eta, double x, double e )
{
  double x_pow_3_2 = sqrt(x) * x;
  double phidot = 0;

  if (conservative_pn_order == 0) /* 0pN term */
  {
    phidot = phi_dot_0pn(e, eta, u) * x_pow_3_2;
  }
  else if (conservative_pn_order == 1) /* 1 pN term */
  {
    phidot = ( phi_dot_0pn(e, eta, u) + x * phi_dot_1pn(e, eta, u) ) * x_pow_3_2;
  }
  else if (conservative_pn_order == 2) /* 2 pN term */
  {
    phidot = ( phi_dot_0pn(e, eta, u) + x * phi_dot_1pn(e, eta, u)
      + x * x * phi_dot_2pn(e, eta, u) ) * x_pow_3_2;
  }
  else if (conservative_pn_order == 3) /* 3 pN term */
  {
    phidot = ( phi_dot_0pn(e, eta, u) + x * phi_dot_1pn(e, eta, u)
      + x * x * phi_dot_2pn(e, eta, u) + 
      x * x * x * phi_dot_3pn(e, eta, u) ) * x_pow_3_2;
  }
  else 
  {
    fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
    exit(1);
  }   

  return phidot;
}
