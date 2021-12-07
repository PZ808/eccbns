/* $Id: relative_separation_guide.c,v 1.1 2009/02/14 16:12:40 pzimmerman Exp $ */

/*
 * 	relative_separation_guide.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>

#define M_SUN 1.98892e30
#define T_SUN 4.92549095e-6
#define L_SUN 1476.62504
#define G_SI 6.67259e-11
#define C_SI 299792458.0
#define DAY 86400.0

double relative_separation( double u );
double kepler_guide( double e, double l, double eps_root );

const double AU = 1.49597870691e11;

/* mass of sun in solar units */
double mass_sun = 1.0;
/* mass of mercury in solar units */
double mass_merc = 1.7634e-7;
double eccentricity = 0.2056;
/* semi-major axis in SI units */
double semi_major_axis = 5.791e10;     
/* period of mercury in SI units */
double period = 7.6005438e6;

double relative_separation( double ecc_anomaly )
{
  return semi_major_axis * ( 1.0 - eccentricity * cos( ecc_anomaly ) );
}

double x_cartesian( double ecc_anomaly )
{
	return semi_major_axis * cos(ecc_anomaly);    
}

double y_cartesian( double ecc_anomaly )
{
	double semi_minor_axis = semi_major_axis * sqrt(1.0 - eccentricity*eccentricity);
	return semi_minor_axis * sin(ecc_anomaly);
}

int main( int argc, char *argv[] )
{
	unsigned long int i;
	unsigned long int i_max;
	unsigned long int out_interval;
	double x_coordinate, y_coordinate;
	double x_au, y_au;
	double total_mass;
	double mean_anomaly, mean_motion, ecc_anomaly;
	double tolerance;
	double r;	
	double l_init;
	double t, dt, t_min, t_max;
	double t_days, out_days;
	double sampling_rate, sampling_interval;

	total_mass = mass_sun + mass_merc;

	if ( argc != 7 )
  {
    fprintf( stderr, "ERROR: incorrect number of arguments\n" );
    fprintf( stderr, "usage: %s t-min t-max initial-mean-anomaly tolerance sampling-rate output-interval\n", argv[0] );
    return 1;
  }

	t_min 			  = atof( argv[1] );
	t_max 			  = atof( argv[2] );
	l_init        = atof( argv[3] );
	tolerance 	  = atof( argv[4] );
	sampling_rate = atof( argv[5] );
	out_interval  = atoi( argv[6] ); 

	t_days = t_max/DAY;
	out_days = out_interval/DAY;
	mean_motion = (2.0 * M_PI * T_SUN) / period;

	/* output file parameters */
	FILE *fp = NULL;
	const int fname_length = 256;
	char fname[fname_length];

  /* create the the file name and file pointer */
	snprintf(fname, fname_length * sizeof(char),
			"merc_newton_%.2f_%2.1e_%.1f_%.2e.txt", 
			t_days, tolerance, sampling_rate, (float)out_interval);
	fp = fopen (fname, "w");

	/* sampling rate [Hz] */
	sampling_interval = 1.0 / sampling_rate;
	/* scale the steps by t_sun */
	dt = sampling_interval/T_SUN;

	t = t_min;
  i_max = (unsigned long int)(t_max);
	t_min /= T_SUN;
	t_max /= T_SUN;

	//fprintf( stderr, "n = %f\n", mean_motion);

	for (i = 1; i <= i_max; i++)
	{
		t = dt * (double)i;
		mean_anomaly = l_init + mean_motion * t;
		ecc_anomaly = kepler_guide( eccentricity, 
				mean_anomaly, tolerance );
  	r = relative_separation(ecc_anomaly);
		x_coordinate = x_cartesian(ecc_anomaly);
		y_coordinate = y_cartesian(ecc_anomaly);

    if ( i % out_interval == 0 ) 
    {
    	x_au = x_coordinate/AU; 
    	y_au = y_coordinate/AU; 
    	r /= AU;
    	t *= T_SUN;

			fprintf( fp, "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e\n", 
					t/DAY, r, x_au, y_au, ecc_anomaly, mean_anomaly ); 
		}
	}	

	fclose( fp );
	return 0;
}
