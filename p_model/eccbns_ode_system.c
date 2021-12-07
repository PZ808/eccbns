/* 
 *  eccbns_ode_system.c							                 
 *  $Id: eccbns_ode_system.c,v 1.3 2009/05/31 18:55:09 pzimmerman Exp $
 * 											 
 *  Program sets up the framework for the ode system
 *  corresponding to the eccentric binary initial value problem. 
 *  											 	
 */ 

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include "ode_system.h"

extern double e_init;
extern double f_init;
extern double m1;
extern double m2;

/* Euler angles (degree measure) */
const double iota_deg  = 0.0;     /* source's first polar angle       */
const double beta_deg  = 0.0;     /* source's second polar angle   	  */
const double theta_deg = 0.0;     /* detector's first polar angle   	*/
const double phi_deg   = 0.0;     /* detector's second polar angle  	*/
const double psi_deg   = 0.0;     /* polarization plane's polar angle */

/* mass of the sun in seconds */
const double t_sun_si = 4.92549095e-06;
/* time required for light to travel one parsec in seconds */
const double t_parsec = 1.0292712503e8; 

/* function: init_y 
 * sets the initial conditions 
 */
int init_y( double t, double y[], void *params )
{
	double *ode_parameters = (double*) params;
	double mass1 = ode_parameters[0];
	double mass2 = ode_parameters[1]; 	
	double f_init = ode_parameters[2];
	double e_init = ode_parameters[3];
	double M = mass1 + mass2;  	

  /* hardwire to t = 0 */
	if ( t == 0 )
	{ 
    /* intitial phase */ 
		y[0] = 0.0; 
    /* initial semi-latus rectum */
		y[1] = (1.0 - e_init*e_init) / pow(M_PI * M * t_sun_si * f_init, 2./3.);	  
    /* initial eccentricity */
		y[2] = e_init; 
	}
	else
	{
		fprintf( stderr, "Error: initial conditions are only known for t = 0\n" );
		exit( 1 );
	}
	return 0;
}

/* function: ode_system
 * sets up the ODE system
 * as given in PRD 60 124008 
 */ 
int ode_system( double t, const double y[], double dydt[], void *params )
{
	double *ode_parameters = (double*) params;
	double mass1 = ode_parameters[0];
	double mass2 = ode_parameters[1];
	double cos_phase = cos( y[0] );
	double p = y[1];
	double e = y[2];
	double M = mass1 + mass2;   
	double mu = mass1*mass2 / (mass1 + mass2); 
  double mass_factor = mu/(M*M);

	/* p-model evolution equations 
   * 0pN eccentric binary system */

	dydt[0] = ( 1.0 + e * cos_phase ) * ( 1.0 + e * cos_phase ) / 
		( M * pow( p, 3.0/2.0 ) );

	dydt[1] = ( -64.0 / 5.0 ) * ( mass_factor ) * (pow( 1.0 - e*e, 3./2. ) /
			( p*p*p ) ) * ( 1.0 + 7.0 * e*e / 8.0 );

	dydt[2] = ( -304.0 / 15.0 ) * ( mass_factor ) * (pow( 1.0 - e*e, 3./2.) / 
			( p*p*p*p ) ) * e * ( 1.0 + 121.0 * e*e / 304.0 );

	return GSL_SUCCESS;
}

int ode_term_cond( double t, double y[], void *params )
{

	if ( y[1] <= 6.0 )
		return 1;
	else
		return 0;
}

char* ode_filename( void )
{
	int name_len = 128;
	char *name = malloc( name_len * sizeof(char) );
	snprintf( name, name_len * sizeof(char), "%.1f_%.1f_%.1f_%.2f", m1, m2, f_init, e_init );
	return name;
}

int write_output( FILE *fp, double t, double t_scale, double y[], void *params  )
{
	double *gw_params = (double*) params;
	double mass1, mass2;
	double m_tot, mu;
	double iota, beta, theta, phi, psi;
	double phase, p, e;
	double f_plus, f_cross;
	double h_plus, h_cross;
  double h_plus_dom, h_cross_dom;
	double h, h_dom;

	const double R = 1.0e-9 * t_parsec; 

	iota  = iota_deg * M_PI / 180.0;
	beta  = beta_deg * M_PI / 180.0;
	theta = theta_deg * M_PI / 180.0;
	phi   = phi_deg * M_PI / 180.0;
	psi   = psi_deg * M_PI / 180.0;

	mass1 = gw_params[0];
	mass2 = gw_params[1];

	m_tot = mass1 + mass2;
	mu    = mass1*mass2 / (m_tot); 

	phase = y[0];
	p     = y[1];
	e     = y[2];

	 /* beam factors for the instrument */
	f_plus   = 0.5 * (1. + pow(cos(theta), 2.)) * cos(2.* phi) * cos(2.* psi) 
		- cos(theta) * sin(2.* phi) * sin(2.* psi);
	f_cross  = 0.5 * (1. + pow(cos(theta), 2.)) * cos(2.* phi) * sin(2.* psi) 
		+ cos(theta) * sin(2.* phi) * cos(2.* psi);
	
	/* waveform polarization equations direct from H. Wahlquist 
	 * Gen. Relativ. Gravit. 19, 1101 (1987) */
	h_plus  = - (mu / (p*R)) * ( ( 2.* cos(2.* phase - 2.* beta) + (5.* e/2.) * cos(phase - 2.* beta) 
				+ (e/2.) * cos(3.* phase - 2.* beta) + e*e * cos(2.* beta) ) 	
			* (1. + cos(iota)*cos(iota)) - e * (cos(phase) + e) * sin(iota)*sin(iota) );

	h_cross = - (mu / (p*R)) * ( 4.* sin(2.* phase - 2.* beta) + 5.* e * sin(phase - 2.* beta)
			+ e * sin(3.* phase - 2.* beta) - 2.* e*e * sin(2.* beta) ) * cos(iota);

	h = f_plus * h_plus + f_cross * h_cross;

  h_plus_dom = -( (2.0*mu) / (p*R) ) * ( (1.0 + cos(iota)*cos(iota)) * 
    cos(2.0*phase-2.0*beta) - e*e * sin(iota)*sin(iota) );
  
  h_cross_dom = - ( (4.0*mu*cos(iota)) / (p*R) ) * (
      sin(2.0*phase-2.0*beta) - 2.0 * e*e * sin(2.0*beta) );

  /* dominant harmonic only */
	h_dom = f_plus * h_plus_dom + f_cross * h_cross_dom;

	fprintf(fp, "%24.16e %24.16e\n",
		  t * t_scale, h_dom, h );
#if 0
	fprintf(fp, "%24.16e %24.16e %24.16e %24.16e %24.16e "
    "%24.16e %24.16e %24.16e %24.16e %24.16e %24.16e %24.16e \n",
		  t * t_scale, t, phase, p, e, h_plus, h_cross, h, 
      amp_plus, amp_cross, h_plus_dom, h_cross_dom );
#endif
	return 0;
}

/* function: check_args */
int check_args( void )
{
	if ( e_init < 0 || e_init > 1) 
	{ 
		fprintf( stderr, "Error: invalid argument to --e-init (%f)\n", e_init ); 
		return 1;
	} 

	if ( f_init <= 0 ) 
	{ 
		fprintf( stderr, "Error: invalid argument to --f-init (%f)\n", f_init ); 
		return 1;
	} 

	if ( m1 < 0 ) 
	{ 
		fprintf( stderr, "Error: invalid argument to --mass-1 (%f)\n", m1 ); 
		return 1;
	} 

	if ( m2 < 0 ) 
	{ 
		fprintf( stderr, "Error: invalid argument to --mass-2 (%f)\n", m2 ); 
		return 1;
	} 
	return 0;
}
