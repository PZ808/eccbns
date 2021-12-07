/* $Id: eT_35pn.c,v 1.4 2009/02/21 18:45:53 pzimmerman Exp $ */

/*
 * eT_35pn.c
 *
 * 	Program numerically evolves a system of ODE's (using GSL integration routines) 
 * 	and generates a 3.5-pN Taylor-Et waveform. 
 *  
 */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h> 
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "eT_ICs.h"


int et_ode_system (double t, const double y[], double dydt[],
		void *params)
{

	double *et_params = (double *) params;
	double mass1 = et_params[0];
	double mass2 = et_params[1];
	//double phi = y[0];
	double xi  = y[1];
	double M = mass1 + mass2;
	double mu = mass1*mass2 / ( mass1 + mass2 );
	double eta = mu /M;
	const double euler_gamma = 0.57721566490153286;
	//const double t_sun = 4.92549095e-6;   /* mass of the sun in seconds   */

	dydt[0] = (pow(xi, 3./2.) / M) * ( 1. + (1./8.) * ( 9. + eta ) * xi
			+ ( 891./128. - (201./64.)* eta + (11./128.) * eta*eta ) * xi*xi
			+ ( 41445./1024. - ( 309715./3072. - (205./64.)* M_PI*M_PI ) * eta
				+ (1215./1024.)* eta*eta + (45./1024.)* pow(eta,3.) ) * pow(xi,3.) );

	dydt[1]	= 64./5.* eta * pow(xi, 5.) * ( 1. + ( 13./336. - (5./2.)* eta ) * xi
			+ 4.* M_PI * pow(xi, 3./2.) + ( 117857./18144. - (12017./2016.)* eta 
				+ (5./2.)* eta*eta ) * xi*xi + ( 4913./672. - (177./8.)*eta ) * M_PI * pow(xi, 5./2.)
			+ ( 37999588601./279417600. - (1712./105.)* log( 4.* sqrt(xi) )
				- (1712./105.)* euler_gamma + (16.* M_PI*M_PI) / 3. 
				+ ( (369./32.)* M_PI*M_PI - 24861497./72576. ) * eta + (488849./16128.)* eta*eta 
				- 85./64.* pow(eta,3.) ) * pow(xi,3.) + ( 129817./2304. - 3207739./48384.* eta 
				+ 613373./12096.* eta*eta) * M_PI * pow(xi, 7./2.) );

	//fprintf( stderr, "dphi_dt = %e\t dxi_dt = %e\n", dydt[0], dydt[1] );

	return GSL_SUCCESS;
}

int main( int argc, char* argv[] )
{
	double mass1, mass2;
	double f_gw;
	double omega;
	double M, mu, eta;
	//double phi;
	//double xi;
	double h;
	//double x, E;
	double N_cycles;

	FILE* fp = NULL;
	FILE* fp_dat = NULL;
	const int fname_len = 256;
	char fname[fname_len];

	/* time stepping variables for integrating ode's */
	int i;
	int status; 
	double t, t_next;
	double step, step_size;
	double samp_rate = 16384.0;
	double dt = 1.0 / samp_rate;

	/* variables for integration */
	double y[2];
	double eT_params[3];
	const double c_si = 299792458.0;        /* m s^{-1}                     */
	const double t_sun = 4.92549095e-6;     /* mass of the sun in seconds   */
	//const double g_si = 6.67259e-11;

	/* gsl solver mechanism */

	const gsl_odeiv_step_type* solver_type
		= gsl_odeiv_step_rkf45;
	gsl_odeiv_step* solver_step
		= gsl_odeiv_step_alloc( solver_type, 2 );
	gsl_odeiv_control* solver_control
		= gsl_odeiv_control_standard_new( 0.0, 1.0e-10, 1.0, 1.0 );
	gsl_odeiv_evolve* solver_evolve
		= gsl_odeiv_evolve_alloc ( 2 );
	gsl_odeiv_system solver_system = { et_ode_system, 
		NULL, 2, 
		(void*) eT_params };

	if ( argc != 4 )
	{
		fprintf( stderr, "error: incorrect number of arguments\n" );
		fprintf( stderr, "usage: %s mass1 mass1 f_gw\n", argv[0] );
		return 1; 
	}

	mass1 = eT_params[0] = atof( argv[1] ); /* mass1 [solar masses] */
	mass2 = eT_params[1] = atof( argv[2] ); /* mass2 [solar masses] */        
	f_gw  = eT_params[2] = atof( argv[3] ); /* initial GW freq. [Hz] */

	M = mass1 + mass2;
	mu = (mass1*mass2) / ( mass1 + mass2);
	eta = mu / M;
	omega = M_PI * f_gw;

	double init_params[3];

	init_params[0] = M;
	init_params[1] = eta;
	init_params[2] = f_gw;     /* initial GW freq. [Hz] */

	//energy_params.M = M;
	//energy_params.eta = eta;


	snprintf( fname, fname_len * sizeof(char), "eT_waveforms_%4.2f_%4.2f.txt", mass1, mass2 );
	fp = fopen( fname, "w" );

	snprintf( fname, fname_len * sizeof(char),
			"eT_waveforms_%4.2f_%4.2f.dat", mass1, mass2 );
	fp_dat = fopen( fname, "w" );
	fprintf( fp_dat, "%% dx = %24.20e\n", dt );

	/*
	 * make the root solving mechanism
	 */

	int root_status;
	int iter = 0;
	int max_iter = 1000;
	const gsl_root_fsolver_type *root_finder_type;
	gsl_root_fsolver *root_solver;
	double xi_root = 0;
	double xi_lo   = 0.; 	         /* estimated lower bracket */
	double xi_hi   = 200.;  	     /* estimated upper bracket */
	double root_eps_abs = 0.0;
	double root_eps_rel = 1.0e-6;

	gsl_function E;

	//struct energy_params rootget_params = { M, eta };

	E.function = &get_xi_init;
	E.params   = &init_params;

	root_finder_type = gsl_root_fsolver_brent;
	root_solver = gsl_root_fsolver_alloc (root_finder_type);

	gsl_root_fsolver_set( root_solver, &E, xi_lo, xi_hi );

	iter = 0;

	do  {
		root_status = gsl_root_fsolver_iterate (root_solver);
		xi_root = gsl_root_fsolver_root (root_solver);
		xi_lo = gsl_root_fsolver_x_lower (root_solver);
		xi_hi = gsl_root_fsolver_x_upper (root_solver);
		root_status = gsl_root_test_interval ( xi_lo, xi_hi, 
				root_eps_abs,
				root_eps_rel);

		if ( root_status == GSL_SUCCESS )
		{	
			fprintf( stderr, "xi_initial = %e\n", xi_root );
			break;
		}

		if ( iter == max_iter || root_status != GSL_CONTINUE )
		{
			fprintf( stderr, "Failure in GSL root finder\n" );
			exit( 1 );
		}
	}  while ( 1 );

	/*
	 * set initial values
	 */

	t    = 0.0;
	y[0] = 0.0;
	y[1] = xi_root;

	fprintf( stderr, "phi_init = %f \t xi_init = %f \n", y[0], y[1] );
	i = 0;

	/* scale the step size */
	step = dt/t_sun;

	while ( 1 )
	{
		/* advance the time */
		t =  i * step;
		t_next = ++i * step;  
		step_size = step;

		status = gsl_odeiv_evolve_apply(
				solver_evolve, solver_control, solver_step, &solver_system,
				&t, t_next, &step_size, y );

		if ( status != GSL_SUCCESS )
		{
			fprintf( stderr, "failure in GSL integrator\n" );
			break;
		}

		/* end evolution when the NANs appear */
		if ( isnan( y[0] ) || isnan( y[1] ) )
		{
			fprintf( stderr, "y[0] is %f, y[1] is %f\n", y[0], y[1] );
			exit( 1 );
		}

		/* generate eT waveform */
		h = - ( 0.5 * mu * y[1] * c_si * c_si  ) * cos( 2.* y[0] );

		fprintf( stderr, "%32.16e %32.16e %32.16e %32.16e\n", t*t_sun, y[0], y[0]/(2.*M_PI), h );
		fprintf( fp, "%32.16e %32.16e %32.16e %32.16e\n", t*t_sun, y[0], y[0]/(2.*M_PI), h );
		fprintf( fp_dat, "%32.16e\n", h );

	}

	fclose( fp );
	fclose( fp_dat );

	gsl_root_fsolver_free (root_solver);
	gsl_odeiv_evolve_free (solver_evolve);
	gsl_odeiv_control_free (solver_control);
	gsl_odeiv_step_free (solver_step);

	return 0;
}
