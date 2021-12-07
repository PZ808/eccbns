/* $Id: x_ecc_kepler.c,v 1.8 2009/05/17 06:17:25 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>

#include "x_ecc_wave.h"

static double cosu_factor(double e, double u)
{
	return (e * cos(u) - 1.);
}

static double ecc_phi_0pn( double e, double eta ) /* Eqn (A21) */
{
	return e;
}

static double ecc_phi_1pn( double e, double eta ) /* Eqn (A22) */
{
	return e * ( 4.0 - eta );
}

static double ecc_phi_2pn( double e, double eta ) /* Eqn (A23) */
{
	double e_pow_2 = e*e;
	return  ( e / (96.* (e_pow_2 - 1.)) ) * 
    ( ( (41.* eta - 659.)*eta + 1152. ) * e_pow_2 +
			(4.* eta + 68.)*eta + sqrt(1.-e_pow_2) * (288.* eta - 720.) - 1248.);
}

static double ecc_phi_3pn( double e, double eta ) /* Eqn (A24) */
{
	double eta_pow_2 = eta*eta;
	double e_pow_2 = e*e;
	double pi_pow_2 = M_PI*M_PI;

	return ( - (e / (26880.* pow(1. - e_pow_2, 5./2.))) *
			( ( (13440.* eta_pow_2 + 483840.* eta - 940800.) * e_pow_2 +
					(255360.* e_pow_2 + (17220.* pi_pow_2 - 2880640.) * eta
					 + 2688000) ) * e_pow_2 +
				(-268800.* eta + 2396800.) * eta
				+ sqrt(1. - e_pow_2) *
				( ( ( (1050.* eta - 134050.)*eta_pow_2 + 786310.* eta - 860160.)*e*e
						+ ((-18900.* eta + 553980.)*eta_pow_2 + 4305.* pi_pow_2*eta
							- 1246368.* eta + 2042880.) ) * e_pow_2 + (276640.*eta +
							2674480. - 17220.* pi_pow_2) * eta - 1451520.)
				- 17220.*M_PI*M_PI*eta - 1747200.) );
}

static double l_0pn( double e, double u, double eta ) /* Eqn (A16) */
{
	return u - e*sin(u);
}

static double l_2pn( double e, double u, double eta,
		double e_phi, double v_minus_u ) /* Eqn (A17) */
{
	double u_factor = cosu_factor(e, u);
	double e_pow_2 = e*e;
	double e_factor = 1.0 - e_pow_2;

	return ( ( 1. / ( 8.* sqrt(e_factor) * (1.0 - e * cos(u)) ) ) *
			( -12.* (2.* eta - 5.) * (-v_minus_u) * u_factor -
				e * sqrt(e_factor) * (eta - 15.) * eta * sin(u) ));
}

static double l_3pn( double e, double u,
		double eta, double e_phi, double v_minus_u ) /* Eqn (A18) */
{
	double u_factor = cosu_factor(e, u);
	double u_minus_v = -v_minus_u;
	double u_factor_pow_3 = u_factor*u_factor*u_factor;
  double neg_u_factor_pow_3 = (1.0-e*cos(u))*(1.0-e*cos(u))*(1.0-e*cos(u));
	double e_pow_2 = e*e;
	double e_factor = (1.0 - e_pow_2);
	double eta_pow_2 = eta*eta;

	return ( ( 1. / (6720.* e_factor * sqrt(e_factor) * neg_u_factor_pow_3) ) *
		( 35.* ( 96.* (11.*eta_pow_2 - 29.* eta + 30.) * e_pow_2
		 	+ 960.* eta_pow_2 + eta * (-13184. + 123.* M_PI*M_PI) + 8640. ) * u_minus_v * (u_factor_pow_3) +
		 	 3360.* ( -12.* (2.* eta - 5.) * (u_minus_v) + 12.* e * (2.* eta - 5.) * cos(u)*u_minus_v + 
          e * sqrt(e_factor) * (eta - 15.) * eta * sin(u) ) * u_factor*u_factor +
		        e * sqrt(e_factor) * (140.* ((13.* e_pow_2 - 11.) * e_pow_2 - 2.) * eta*eta_pow_2 -
					    140.* ( (73.* e_pow_2 - 325.) * e_pow_2 + 444.) * eta_pow_2 +
					     ( (3220.* e_pow_2 - 148960.) * e_pow_2 - 4305.* M_PI*M_PI + 143868.) * eta +
					 e_pow_2 * ( 1820.* (e_pow_2 - 1.) * eta*eta_pow_2 - 140. * (83.* e_pow_2 + 109.) * eta_pow_2 - 
             ( 1120.* e_pow_2 + 4305.* M_PI*M_PI + 752. ) * eta + 67200. ) * cos(u)*cos(u) - 
           2.* e * ( 1960.* (e_pow_2 - 1.) * eta*eta_pow_2 + 6720.* (e_pow_2 - 5.) * eta_pow_2 + 
             ( -71820.* e_pow_2 - 4305.* M_PI*M_PI + 69948. ) * eta + 67200. ) * cos(u) + 67200.) * sin(u) ) );
}

static double ecc_phi( int conservative_pn_order,
		double eta, double x, double e )  /* Eqn. (A21) */
{
	double ephi = ecc_phi_0pn( e, eta );

	if ( (conservative_pn_order < 0) || (conservative_pn_order > 3)  )
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit (1);
	}
	if ( conservative_pn_order > 0 )      /* 1pN term */
		ephi += ecc_phi_1pn( e, eta ) * x;
	if ( conservative_pn_order > 1 )      /* 2pN term */
		ephi += ecc_phi_2pn( e, eta ) * x*x;
	if ( conservative_pn_order > 2 )      /* 3pN term */
		ephi += ecc_phi_3pn( e, eta ) * x*x*x ;

	return ephi;
}

static double mean_anomaly( int conservative_pn_order,
		double eta, double u, double x, double e )  /* Eqn. (A25) */
{
	double x_pow_2 = x*x;
	double e_phi = ecc_phi( conservative_pn_order, eta, x, e );  /* Eqn. (A21) */
	double beta_phi =  ( (1.0 - sqrt(1.0 - e_phi*e_phi)) / e_phi );
	double v_minus_u = ( 2.0 * atan( (sin(u) * beta_phi)) /
			(1.0 - beta_phi * cos(u)) ) ;

	double l = l_0pn( e, u, eta );

	if ( (conservative_pn_order < 0) || (conservative_pn_order > 3)  )
	{
		fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
		exit (1);
	}

	if ( conservative_pn_order > 1 ) /* 2pN term */
		l += l_2pn( e, u, eta, e_phi, v_minus_u ) * x_pow_2;
	if ( conservative_pn_order > 2 ) /* 3pN term */
		l += l_3pn( e, u, eta, e_phi, v_minus_u ) * x * x_pow_2;

	return l;
}

double kepler_equation( double u, void* params )
{
	struct kepler_params *k = (struct kepler_params *) params;
	int pn_order = k->pn_order;
	double result;

	if ( pn_order < 2 )
	{
		result = k->l - (u - (k->e) * sin(u));
	}
	else
	{
		result = k->l - mean_anomaly( pn_order, k->eta, u, k->x, k->e);
	}

	return result;
}

double pn_kepler_equation( int conservative_pn_order, double eta, double x,
		double e, double l )
{
	int kep_root_status = 0;
	double u = 0.0;
	double u_expected = mikkola_finder( e, l );
	unsigned int iter = 0;
	unsigned int max_iter = 100;
	double u_low = M_PI;  /* left bracket */
	double u_high = M_PI;  /* right bracket */
	double root_eps_abs = 0.0;
	double root_eps_rel = 1.0e-3;

  double sgn_l;
	const gsl_root_fsolver_type *root_finder_type;

	gsl_root_fsolver *root_solver;
	gsl_function kepler_gsl_func;

	root_finder_type = gsl_root_fsolver_brent;
	root_solver = gsl_root_fsolver_alloc( root_finder_type );

	/* zero mean anomaly case */
	if ( (conservative_pn_order < 0) || (conservative_pn_order > 3) )
	{
		exit (1);
		fprintf( stderr, "error: error in pN order\n" ); 
	}

	if ( l == 0 )
	{	
		return u;
	}

	if ( conservative_pn_order < 2 ) 
  {
		return u_expected;
  }

	if ( conservative_pn_order == 2 ) 
	{

    /* range reduction of mean_anomaly  */
    while ( l > M_PI )
    {
      l -= 2*M_PI;
    }
    while ( l < -M_PI )
    {
      l += 2*M_PI;
    }

    /* compute the sign of l */
    if ( l >= 0.0 )
      l = 1.0;
    else
    {
      l = -1.0;
    }

    l *= sgn_l;

	  struct kepler_params params = {conservative_pn_order, eta, x, e, l};
	  kepler_gsl_func.function = &kepler_equation;
	  kepler_gsl_func.params = (void*) &params;
    gsl_root_fsolver_set( root_solver, &kepler_gsl_func, u_low, u_high );

		do {

			iter++;

			kep_root_status = gsl_root_fsolver_iterate( root_solver );
			u = gsl_root_fsolver_root( root_solver );
      fprintf( stderr, "u=%2.2f\n",u);
      u_low = gsl_root_fsolver_x_lower(root_solver);
      u_high = gsl_root_fsolver_x_upper(root_solver);
			kep_root_status = gsl_root_test_interval( u_low, u_high,
					root_eps_abs, root_eps_rel );

			if ( kep_root_status == GSL_SUCCESS ) 
      {
	      u *= sgn_l;
	      return u;
      }
      if ( iter == max_iter || kep_root_status != GSL_CONTINUE )
			{
				fprintf( stderr, "Error in GSL root finder: %d \n", kep_root_status );
				exit( 1 );
			}
		} while ( 1 );
	}

	gsl_root_fsolver_free( root_solver );

	return u;
}

/* function to solve Keppler's equation and return the eccentric anomaly */
/* inputs are the eccentricity and the mean anomaly                      */
double mikkola_finder( double eccentricity, double mean_anomaly )
{
	/* variables needed for Mikkola's method */
	double a, b, sgn_b, z, s;
	double sgn_mean_anomaly;
	double ecc_anomaly;

	/* range reduction of mean_anomaly  */
	while ( mean_anomaly > M_PI )
	{
		mean_anomaly -= 2*M_PI;
	}
	while ( mean_anomaly < -M_PI )
	{
		mean_anomaly += 2*M_PI;
	}

	/* compute the sign of l */
	if ( mean_anomaly >= 0.0 )
		sgn_mean_anomaly = 1.0;
	else
	{
		sgn_mean_anomaly = -1.0;
	}

	mean_anomaly *= sgn_mean_anomaly;

	/* compute alpha and beta of Minkkola Eq. (9a) */
	a  = ( 1.0 - eccentricity ) / ( 4.0 * eccentricity + 0.5 );
	b  = ( 0.5 * mean_anomaly ) / ( 4.0 * eccentricity + 0.5 );

	/* compute the sign of beta needed in Eq. (9b) */
	if ( b >= 0.0 )
		sgn_b = 1.0;
	else
		sgn_b = -1.0;

	/* Mikkola Eq. (9b) */
	z = pow( ( b + sgn_b * sqrt(b*b + a*a*a) ), 1./3. );
	/* Mikkola Eq. (9c) */
	s = z - a / z;
	/* add the correction given in Mikkola Eq. (7) */
	s = s - 0.078 * pow(s, 5.0) / ( 1.0 + eccentricity );
	/* finally Mikkola Eq. (8) gives u */
	ecc_anomaly = mean_anomaly + eccentricity *
		( 3.0 * s - 4.0 * pow(s, 3.0) );
	/* correct the sign of u */
	ecc_anomaly *= sgn_mean_anomaly;

	return ( ecc_anomaly );
}
