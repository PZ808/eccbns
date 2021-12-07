/* $Id: mikkola.c,v 1.2 2009/02/18 02:17:33 pzimmerman Exp $ */

#include <math.h>
#include <gsl/gsl_math.h>

double mikkola_finder( double e, double l); 

double mikkola_finder( double eccentricity, 
		double mean_anomaly )
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

