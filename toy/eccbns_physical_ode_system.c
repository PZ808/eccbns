/* 
 *  eccbns_physical_ode_system.c							                 
 * 											 
 *  Program sets up the framework for the ode system
 *  corresponding to the eccentric 0pN binary initial 
 *  value problem. Runs our code in "Eric Mode".  
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

/* Euler angles: degrees */

const double iota_deg  = 45.0;      /* source's first polar angle    	*/
const double beta_deg  = 0.0;      /* source's second polar angle   	*/
const double theta_deg = 45.0;      /* detector's first polar angle   	*/
const double phi_deg   = 20.0;      /* detector's second polar angle  	*/
const double psi_deg   = 0.0;       /* polarization plane's polar angle */


const double t_parsec = 1.0292712503e8; 



int init_y( double t, double y[], void *params )
{
	double *ode_parameters = (double*) params;
	double mass1 = ode_parameters[0];
	double mass2 = ode_parameters[1]; 	
	double f_init = 26.6666;
	double e_init = ode_parameters[3];
	double M = (mass1 + mass2)*4.92549095e-06;  	

	if ( t == 0 )
	{
		y[0] = 0.0;
		y[1] = ( 1. - e_init*e_init ) / pow( M_PI * M * f_init, 2./3. );
		y[2] = e_init;

	}
	else
	{
		fprintf( stderr, "Error: initial conditions are only known for t = 0\n" );
		exit( 1 );
	}

	fprintf( stderr, "M_PI   = %32.16e\n", M_PI );
	fprintf( stderr, "M      = %32.16e\n", M );
	fprintf( stderr, "f_init = %32.16e\n", f_init );
	fprintf( stderr, "p_init = %32.16e\n", y[1] );
	fprintf( stderr, "e_init = %32.16e\n", y[2] );

	return 0;
}

int ode_system( double t, 
		const double y[], 
		double dydt[],
		void *params )
{
	double *ode_parameters = (double*) params;
	double mass1 = ode_parameters[0];
	double mass2 = ode_parameters[1];
	double cos_phase = cos( y[0] );
	double p = y[1];
	double e = y[2];
	double M = (mass1 + mass2)*4.92549095e-06;            /* total mass in units of solar mass */ 
	double mu = (mass1*mass2 / (mass1 + mass2))*4.92549095e-06; 

	/* evolution equations for a 0 pN binary system */

	dydt[0] = ( 1.0 + e * cos_phase ) * ( 1.0 + e * cos_phase ) / 
		( M * pow( p, 3.0/2.0 ) );

	dydt[1] = ( -64.0 / 5.0 ) * ( mu / (M*M) ) * (pow( 1.0 - e*e, 3.0 / 2.0 ) /
			( p*p*p ) ) * ( 1.0 + 7.0 * e*e / 8.0 );

	dydt[2] = ( -304.0 / 15.0 ) * ( mu / (M*M) ) * (pow( 1.0 - e*e, 3.0/2.0) / 
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
	snprintf( name, name_len * sizeof(char), "%.2f_%.2f_%.2f", m1, m2, e_init );
	return name;
}

int write_output( FILE *fp, double t, double t_scale, double y[], void *params  )
{
	double *gw_params = (double*) params;
	double mass1;
	double mass2;
	double mu;
	double iota, beta, theta, phi, psi;
	double R    = 1.0e-9 * t_parsec; 
	double phase;
	double p;
	double e;
	double f_plus, f_cross;
	double h_plus, h_cross;
	double h;

	iota  = iota_deg * M_PI / 180.0;
	beta  = beta_deg * M_PI / 180.0;
	theta = theta_deg * M_PI / 180.0;
	phi   = phi_deg * M_PI / 180.0;
	psi   = psi_deg * M_PI / 180.0;

	mass1 = gw_params[0];
	mass2 = gw_params[1];
	mu    = mass1*mass2 / (mass1 + mass2); 

	phase = y[0];
	p   = y[1];
	e   = y[2];

	f_plus   = 0.5 * (1. + pow( cos(theta), 2.)) * cos(2. * phi) * cos(2. * psi) 
		- cos(theta) * sin( 2. * phi ) * sin( 2. * psi );

	f_cross  = 0.5 * ( 1. + pow( cos(theta), 2.) ) * cos( 2. * phi ) * sin( 2. * psi ) 
		+ cos(theta) * sin( 2. * phi ) * cos( 2. * psi );
	
	/*
	 *
	 * These are the correct polarization equations direct from H. Wahlquist
	 * Gen. Relativ. Gravit. 19, 1101 (1987) 
	 *
	 */

	h_plus  = - (mu /(p*R)) * ( ( 2.* cos(2.* phase - 2.* beta) + (5.* e/2.) * cos(phase - 2.* beta) 
				+ (e/2.) * cos(3.* phase - 2.* beta) + e*e * cos(2.* beta) ) 	
			* (1. + cos(iota)*cos(iota)) - (e*cos(phase) + e*e) * sin(iota)*sin(iota) );

	h_cross = - (mu/(p*R)) * ( 4.* sin(2.* phase - 2.* beta) + 5.* e * sin(phase - 2.* beta)
			+ e * sin(3.* phase - 2.* beta) - 2.* e*e * sin(2.* beta) ) * cos(iota);    

	h = f_plus * h_plus + f_cross * h_cross;

	fprintf( fp, "%32.16e %32.16e %32.16e %32.16e %32.16e %32.16e %32.16e %32.16e \n", t * t_scale, t, phase, p, e, h_plus, h_cross, h );


	/* 
	 *
	 * The equations below are those used by Eric in wave.c 
	 * They were also implemented in the original eccbns_ode_system.c 
	 * 
	 * 			*** THEY ARE INCORRECT ***
	 */        

#if 0	
	h_plus =  (mu / (p * R)) * ( ( 2. * cos(2. * phase - 2.* beta) + (5.* e / 2.) * cos(phase - 2.* beta) 
				+ (e / 2.) * cos(3.* phase - 2. * beta) + e*e * cos(2. * beta) ) 	
			* ( 1. + cos(iota)*cos(iota) ) + (e * cos(phase) + e*e ) * sin(iota)*sin(iota) );
	h_cross  =  ( mu / ( p * R ) ) * ( 4. * sin( 2. * phase - 2. * beta ) + 5. * e  * sin( phase  - 2.* beta ) + 
			e * sin( 3. * phase - 2. * beta ) - 2. * e*e * sin( 2. * beta ) ) * cos(iota);
#endif

	return 0;
}

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



