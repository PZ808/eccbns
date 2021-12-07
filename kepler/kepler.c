#include <math.h>
#include <gsl/gsl_math.h>
#include "kepler.h"
	
double kepler_guide( 
    double eccentricity, 
    double mean_anomaly,
    double tolerance )
{
  int mean_anom_negative = 0;
  double newton_error;
  double newton_threshold;
  int newton_iterations = 0;
  double ecc_anomaly;

  /* zero mean anomaly case */
  if ( mean_anomaly == 0 )
  {
    ecc_anomaly = 0;
    return ecc_anomaly;
  }

  /* range reduction of the mean_anomaly */
  while ( mean_anomaly > M_PI )
  {
  	mean_anomaly -= 2*M_PI;
  }
  while ( mean_anomaly < -M_PI )
  {
  	mean_anomaly += 2*M_PI;
  }

  /* solve in the positive part of the orbit */
  if ( mean_anomaly < 0.0 )
  {
    mean_anomaly = -mean_anomaly;
    mean_anom_negative = 1;
  }


  newton_threshold = tolerance * fabs( 1.0 - eccentricity );
  ecc_anomaly = mean_anomaly;

  /* high eccentricty case */
  if ( (eccentricity > 0.8) && (mean_anomaly < M_PI/3.0) )
  {
    double trial = mean_anomaly / fabs( 1.0 - eccentricity );
    if ( trial*trial > 6.0 * fabs( 1.0 - eccentricity ) )
    {
      /* cubic term is dominant */
      if ( mean_anomaly < M_PI )
        trial = pow( 6.0 * mean_anomaly, 1.0/3.0 );
    }
    ecc_anomaly = trial;
  }

  /* iterarte using Newton's method to get solution */
  if ( eccentricity < 1.0 )
  {
  	newton_error = ecc_anomaly - 
    	eccentricity * sin( ecc_anomaly ) - mean_anomaly;
  	while ( fabs( newton_error ) > newton_threshold )
  	{
    	++newton_iterations;
    	ecc_anomaly -= newton_error / 
      	( 1.0 - eccentricity * cos ( ecc_anomaly ));
    	newton_error = ecc_anomaly - 
      	eccentricity * sin ( ecc_anomaly ) - mean_anomaly;
    }  
  }

  return ( mean_anom_negative ? -ecc_anomaly : ecc_anomaly );
}
