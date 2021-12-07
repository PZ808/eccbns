#include <stdio.h>

double kepler_guide( 
	double eccentricity,
	double mean_anomaly,
	double tolarance);

double relative_separation( 
		double eccentric_anomaly);

double relative_dphi_dt(
		double mean_motion,
		double eccentric_anomaly);

double x_cartesian( 
		double eccentric_anomaly );

double y_cartesian( 
		double eccentric_anomaly );

