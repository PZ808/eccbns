/* $Id: conservative_ode_system.c,v 1.2 2009/03/25 22:14:05 pzimmerman Exp $ */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "conservative_dynamics_pure.h"
#include "conservative_ode_system.h"

/* mass scaled time derivative of the orbital phase */ 
double dphi_dt(
		double ecc_anomaly,
		double eta,
		double e,
		double x )
{
  return (
  ( rel_dphi_dt_0pn( e, ecc_anomaly, eta ) +
  	rel_dphi_dt_1pn( e, ecc_anomaly, eta ) * x +
 		rel_dphi_dt_3pn( e, ecc_anomaly, eta ) * x*x*x ) * x * sqrt(x) );
}
  	
#if 0
	double phi_dot_0pn = rel_dphi_dt_0pn(
	  	e, ecc_anomaly, eta );
	double phi_dot_1pn = rel_dphi_dt_1pn(
	    e, ecc_anomaly, eta );
	double phi_dot_2pn = rel_dphi_dt_2pn(
	    e, ecc_anomaly, eta );

	/* mass scaled dphi/dt */
	double M_times_phi_dot = (
	    ( phi_dot_0pn + phi_dot_1pn * x
	      + phi_dot_2pn * x*x ) * x*sqrt(x) );
	return M_times_phi_dot;
}
#endif

/* mass scaled orbtial separation */
double separation(
		double ecc_anomaly,
		double eta,
		double e,
		double x)
{
  return (
  	rel_sep_0pn( e, ecc_anomaly ) +
  	rel_sep_1pn( e, ecc_anomaly, eta ) / x +
  	(	rel_sep_2pn( e, ecc_anomaly, eta ) +
  		rel_sep_3pn( e, ecc_anomaly, eta ) * x ) * x*x );
}

#if 0
	double r_0pn = rel_sep_0pn(
	  	e, ecc_anomaly );
	double r_1pn = rel_sep_1pn(
	  	e, ecc_anomaly, eta );
	double r_2pn = rel_sep_2pn(
	  	e, ecc_anomaly, eta );

	/* mass scaled separation */
	double r_over_M = r_0pn * (1./x) + r_1pn + r_2pn * x;

	return r_over_M;
}
#endif

/* mass scaled angular eccentrcity */
double angular_eccentricity( 
		double ecc_anomaly,
		double eta,
		double e,
		double x )
{
	return ( e + 
		( angular_ecc_1pn( e, eta ) +
	   	angular_ecc_2pn( e, eta ) * x +
 		 	angular_ecc_3pn( e, eta ) * x*x ) * x );
}                                      

/* mass scaled mean_motion */
double mean_motion_func( 
		double e,
		double eta,
		double x )
{
  return ( 
  	( 1.0 + 
  		mean_motion_1pn( e ) * x +
  		mean_motion_2pn( e, eta ) * x*x +
  		mean_motion_3pn( e, eta ) * x*x*x ) * x * sqrt(x) );
}

#if 0
	double n_1pn = mean_motion_1pn( e );
	double n_2pn = mean_motion_2pn( e, eta );
	/* mass scaled mean_motion */
	return ( (1.0 + n_1pn * x + n_2pn * x*x) * x*sqrt(x) );
}
#endif

#if 0

/* mass scaled mean anomaly */
double mean_anomaly_func(
		double ecc_anomaly,
		double eta,
		double e,
		double x )
{
	double e_phi = angular_eccentricity( ecc_anomaly, eta, e, x );

	double l_0pn = ecc_anomaly - e*cos(ecc_anomaly);
	double l_2pn = mean_anomaly_2pn( e, ecc_anomaly,
			eta, e_phi );
	double l_3pn = mean_anomaly_3pn( e, ecc_anomaly,
			eta, e_phi );

	return ( l_0pn + ( l_2pn + l_3pn * x ) * x*x );
}

#endif
