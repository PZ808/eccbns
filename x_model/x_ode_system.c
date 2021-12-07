#include <stdio.h>
#include <math.h>

/* 
 * Conservative PN Differential Equations
 */

/* 2-PN angular frequency dphi/dt */
double dphi_dt(
		int conservative_pn_order,
		double ecc_anomaly,
		double eta,
		double x,
		double eccentricity )
{
	/* mass scaled dphi/dt */
	if (conservative_pn_order == 0)
	  return( rel_dphi_dt_0pn( eccentricity, 
	  			ecc_anomaly, eta ) * x*sqrt(x) );
	else if (conservative_pn_order == 1)
	  return( ( rel_dphi_dt_0pn( eccentricity, 
	  				ecc_anomaly, eta ) +
					rel_dphi_dt_1pn( eccentricity,
	    			ecc_anomaly, eta ) * x ) * x*sqrt(x) );
	else if (conservative_pn_order == 2)
	  return ( ( rel_dphi_dt_0pn( eccentricity,
	  				ecc_anomaly, eta ) +
					rel_dphi_dt_1pn( eccentricity,
	    			ecc_anomaly, eta ) * x 
	 				rel_dphi_dt_2pn( eccentricity,
	 					ecc_anomaly, eta ) * x*x ) * x*sqrt(x) );
	else
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit(1);
	}
}

/* 2-PN separation */
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

double mean_motion(
		int conservative_pn_order,
		double eta,
		double x,
		double eccentricity)
{
	/* mass scaled mean_motion */
	double n_0pn = pow(x, 3./2.);
	double n_1pn = mean_motion_1pn(eccentricity);
	double n_2pn = mean_motion_1pn(eccentricity);

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

int eccentric_x_model_odes( double t, const double y[],
		double dydt[], void *params )
{
	struct ode_parameters *op = (struct ode_parameters *)params;
	double eta = op->eta;
	double ecc_anomaly = op->ecc_anomaly;
	int conservative_pn_order = op->conservative_pn_order;
	int radiation_pn_order = op->radiation_pn_order;

	/* mass scaled radiative ODEs 
	 *    * dxdt = M*x_dot (A25)
	 *       * dedt = M*e_dot (A30)
	 *          */
	double dxdt;
	double dedt;
	/* mass scaled radiative ODEs 
	 *    * dldt = M*l_dot (8)
	 *       * dphidt = M*phi_dot (A10)
	 *          */
	double dldt;
	double dphidt;

	double x   = y[0];
	double e   = y[1];
	double l   = y[2];
	double phi = y[3];

  if ( radiation_pn_order == 0 )
  {
  	dxdt = x_dot_0pn(e,	eta) * x*x*x*x*x;
  	dedt = e_dot_0pn(e,	eta)* x*x*x*x;
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
  dphidt = dphi_dt( conservative_pn_order, ecc_anomaly, eta, x, e);

  dydt[0] = dxdt;
  dydt[1] = dedt;
  dydt[2] = dldt;
  dydt[3] = dphidt;

  return GSL_SUCCESS;
}

