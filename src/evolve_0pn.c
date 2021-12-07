/* code to evolve 0pn eccentric binaries */
/* $Id: evolve_0pn.c,v 1.16 2008/10/03 03:02:44 pzimmerman Exp $ */
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>

/*  Program: evolve_0pn.c 
 *  
 *	Input Parameters: M [solar masses], eta [dimensionless], 
 *  	n_init [Hz], e_init [dimensionless].
 *
 *  Output: t [sec], n(t) [Hz], e(t), phi(t) [radians], h_cross  
 *
 *  Code originally designed to reproduce the 0pN results of             
 *  single precision results of E. Poisson and K. Martel 
 *  Phys. Rev. D60 (1999) 124008)
 * 
 *  This code runs at double precision numerical accuracy, 
 *  the results should be closer to those of Brown and Zimmerman (2008)
 *  
 *  We adopt the Keplerian formalism as in Gopakumar and Tessmer                                  
 *  arxiv:gr-qc/0712.3199v1                    
 */

int eccentric_binary_evolution (double t, const double y[], double f[],
		void* params)
{
	/* evolution equations for 0 pN binary */
	const double t_sun = 4.92549095e-6;     /* mass of the sun in seconds */
	double* mass_params = (double*) params;  
	double M = mass_params[0];
	double eta = mass_params[1];
	double n = y[0];
	double e = y[1];

	/* f[0] = dn/dt */
	f[0] =  ( pow( M * t_sun * n, 5./3.) * n*n * eta * ( 96. + 292. * e*e + 37. * e*e*e*e ) )  / 
		( 5. * pow( ( 1. - e*e ), 7./2. ) );

	/* f[1] = de/dt */
	f[1] = - ( ( pow( M * t_sun * n, 5./3. ) * n * eta * e ) * ( 304. + 121. * e*e ) ) / 
		( 15. * pow( ( 1. - e*e ), 5./2. ) );

	return GSL_SUCCESS;
}

int main ( int argc, char *argv[] )
{


	/*
	 *
	 * input and output variables
	 *
	 */


	/* input binary parameters */
	double M;     /* total mass           */
	double eta;   /* symmetric mass ratio */
	double n;     /* mean motion = 2pi/T  */
	double e;     /* eccentricity         */
	double p;	  /* semi-latus 		  */ 

	/* output file parameters */
	FILE* fp = NULL;
	FILE* fp_dat = NULL;
	const int fname_len = 256;
	char fname[fname_len];

	/* imput parameters to solver control */
	double eps_abs = 0.0;
	double eps_rel = 1.0e-16;
	double a_y = 1.0;
	double a_dydt = 1.0;


	/*
	 *
	 * variables needed to integrate equations of motion
	 *
	 */


	/* time stepping variables for integrating ode's */
	int i;
	double srate = 4096.0;         /* sample rate (Hz)                    */
	double dt = 1.0/srate;         /* sampling interval (s)               */
	double t, t_next;              /* time index variables                */
	double step, step_size;        /* integrator step size and result     */
	double n_final;                /* final mean motion                   */

	/* dynamical equations for integration */
	double y[2];                  /* dynamical variables to be integrated */
	double mass_params[2];        /* parameters of evolution equations    */

	/* conversion to and from dimensionless units */
	const double c_si = 299792458;        /* m s^{-1}                     */
	const double g_si = 6.67259e-11;      /* kg^{-1} m^3 s^{-2}           */
	const double t_sun = 4.92549095e-6;   /* mass of the sun in seconds   */
	const double m_sun = 1.98892e30;      /* kg                           */


	/*
	 *
	 * variables needed to solve kepler's equation 
	 *
	 */


	/* variables needed for Mikkola's method */
	double a, b, sgn_b, z, s;
	/* mean anomaly: l = f(n) + pn corrections */
	double l, sgn_l;
	/* dimmensionless parameter introduced to simplify things */
	double beta;
	/* true anomaly: v - u = 2 * atan( (beta * sin( u ) ) / ( 1 - beta * cos(u) ) )  */
	double v_minus_u;
	/* eccentric anomaly */
	double u;


	/*
	 *
	 * variables needed to generate polarizations of waveform
	 *
	 */


#if 0
	double h_plus;            /* output plus polarizaion              */
#endif
	double h_cross;           /* output cross polarizaion             */
	double r;       		  /* r = r(u(l), n, e)                    */ 
	double dr_dt;             /* rdot = rdot(u(l), n, e)              */
	double phi;     		  /* phi = phi(u(l), n, e)                */
	double dphi_dt; 		  /* phidot = (u(l), n, e)                */


	/* distance normalization */
	double R = 1.0e-30;


	/*
	 *
	 * make the diff eqn solving mechanism
	 *
	 */


	/* Solver step type is Embedded Runge-Kutta-Fehlberg (4,5) */
	const gsl_odeiv_step_type* solver_type 
		= gsl_odeiv_step_rkf45;
	gsl_odeiv_step* solver_step 
		= gsl_odeiv_step_alloc( solver_type, 2 );
	gsl_odeiv_control* solver_control
		= gsl_odeiv_control_standard_new( eps_abs, eps_rel, a_y, a_dydt );
	gsl_odeiv_evolve* solver_evolve
		= gsl_odeiv_evolve_alloc( 2 );
	gsl_odeiv_system solver_system = { eccentric_binary_evolution,
		NULL, 2, (void*) mass_params };


	/* 
	 *
	 * parse the command line arguments and open the output file
	 *
	 */


	if ( argc != 5 )
	{
		fprintf( stderr, "error: incorrect number of arguments\n" );
		fprintf( stderr, "usage: %s M eta n e\n", argv[0] );
		return 1;
	} 

	M   = mass_params[0] = atof( argv[1] ); /* mass [solar masses]          */
	eta = mass_params[1] = atof( argv[2] ); /* mass ratio [dimensionless]   */
	n   = atof( argv[3] );                  /* mean motion [rad/s]          */
	e   = atof( argv[4] );                  /* eccentricity [dimensionless] */

	fprintf( stderr, "M = %4.2f\n", M );
	fprintf( stderr, "eta = %4.2f\n", eta );
	fprintf( stderr, "n = %4.2f\n", n );
	fprintf( stderr, "e = %4.2f\n", e );

	/* create the file name and file pointer */
	snprintf( fname, fname_len * sizeof(char), 
			"Kep_0pn_%4.2f_%4.2f_%4.2f_%4.2f.txt", M, eta, n, e );
	fp = fopen( fname, "w" );
	snprintf( fname, fname_len * sizeof(char), 
			"Kep_0pn_%4.2f_%4.2f_%4.2f_%4.2f.dat", M, eta, n, e );
	fp_dat = fopen( fname, "w" );
	fprintf( fp_dat, "%% dx = %24.20e\n", dt );


	/*
	 *
	 * set initial and final values
	 *
	 */


	t = 0.0;
	i = 0;
	step = dt;
	y[0] = n;
	y[1] = e;
	n_final = 1. / ( pow( 6.0, 3./2. ) * M * t_sun );

#if 0
	fprintf( stderr, "i = %d, t = %f, f_gw = %e, n = %f, e = %f, tstep = %e\n",
			i, t, y[0] / M_PI, y[0], y[1], step );
#endif

	/*
	 *
	 * generate waveform
	 *
	 */


	/* start computing the waveform we use a while() loop  */
	while ( 1 )
	{
		/* advance the time */ 
		t = i * step;   
		t_next = ++i * step;
		step_size = step;

		/* check for not a number in dynamical variables */

		if ( isnan( y[0] ) || isnan( y[1] ) )
		{
			fprintf( stderr, "y[0] is %f, y[1] is %f\n", y[0], y[1] );
			exit( 1 );
		}

		/* call the solver to get the next time step */

		int status = gsl_odeiv_evolve_apply(
				solver_evolve, solver_control, solver_step, &solver_system,
				&t, t_next, &step_size, y );

		if ( status != GSL_SUCCESS )        
		{ 
			fprintf( stderr, "failure in GSL integrator\n" );
			break;
		}


		if ( y[0] > n_final )
		{
			fprintf( stderr, "breaking due to (n = %e) > (n_final = %e)\n", y[0], n_final );
			break;       
		}


		/*
		 *
		 * Use Minkkola's method to solve Kepler's equation
		 *
		 */


		/* compute the mean anomaly from the mean motion */
		l = y[0] * t;
#if 0
		fprintf( stderr, "%d: t = %e, n = %e, l = %e, ", i, t, n, l );
#endif
		/* range reduction of l */
		while ( l > M_PI )
		{
			l -= 2*M_PI;
		}
		while ( l < -M_PI )
		{
			l += 2*M_PI;
		}

#if 0
		fprintf( stderr, "[l] = %e\n", l );
#endif 

		/* compute the sign of l */
		if ( l >= 0.0 )
			sgn_l = 1.0;
		else
		{
			sgn_l = -1.0;
		}

		l *= sgn_l;

		/* compute alpha and beta of Minkkola Eq. (9a) */
		a  = ( 1. - y[1] ) / ( 4. * y[1] + 0.5);
		b  = ( 0.5 * l ) / ( 4. * y[1] + 0.5 );

		/* compute the sign of beta needed in Eq. (9b) */
		if ( b >= 0.0 )
			sgn_b = 1.0;
		else
			sgn_b = -1.0;

		/* Minkkola Eq. (9b) */
		z = pow( ( b + sgn_b * sqrt( b*b + a*a*a ) ), 1./3. );

		/* Minkkola Eq. (9c) */
		s = z - a / z; 

		/* Add the correction given in Mikkola Eq. (7) */
		s = s - 0.078 * pow( s, 5. ) / ( 1. + y[1] );

		/* Finally Mikkola Eq. (8) gives u */
		u = l + y[1] * ( 3. * s - 4.* pow( s, 3. ) ) ;   

		/* correct the sign of u */
		u *= sgn_l;


		/*
		 *
		 * compute the variable needed to get the waveform
		 *
		 */


		beta      = ( 1. - sqrt(1. - y[1]*y[1]) ) / y[1];
		v_minus_u =  2. * atan( ( beta * sin(u) ) / 
				( 1. - beta * cos(u) ) );  

		r       = pow( (g_si * M * m_sun) / (y[0] * y[0]), 1./3.) * 
			( 1. - y[1] * cos(u) );
		dr_dt   = ( pow( g_si * M * m_sun * y[0], 1./3. ) * 
				( y[1] * sin(u) ) )  / ( 1. - y[1] * cos(u) );

		phi     =  y[0] * t + v_minus_u + y[1] * sin(u);
		dphi_dt = ( y[0] * sqrt( 1. - y[1]*y[1] ) ) / 
			pow( ( 1. - y[1] * cos(u) ), 2. ); 
		p  		= (1. - y[1]*[1] ) / pow( y[0] * M, 2./3. );
#if 0
		h_plus = - ( eta * M / R ) * (  ( 1. + cos(iota)*cos(iota) ) 
				* ( ( M / r + r*r * dphi_dt*dphi_dt - dr_dt*dr_dt ) 
					* cos( 2. * phi ) + 2. * r * dr_dt * dphi_dt * sin( 2. * phi) )  
				+ sin(iota)*sin(iota) * ( M / r - r*r * dphi_dt*dphi_dt - dr_dt*dr_dt ) );                  
#endif

		h_cross = - ( ( 2. * g_si * M * eta ) / ( c_si*c_si*c_si*c_si * R ) ) * 
			( ( (g_si * M) / r + r*r * dphi_dt*dphi_dt - dr_dt*dr_dt ) * sin(2. * phi)  
			  - 2. * dr_dt * r * dphi_dt * cos( 2. * phi ) );

#if 0
		fprintf( stderr, "t=%f, f_gw=%e, n=%f, e=%f, a=%e, b=%e, u=%e\n",
				t, y[0] / M_PI, y[0], y[1], a, b, u );
		fprintf( fp, "%e %e %e %e %e\n", t, y[0], y[1], h_plus, h_cross);
#endif
		fprintf( fp, "%24.20e %24.20e %24.20e %24.20e %24.20e %24.20e\n", t, phi, p, y[0], y[1], h_cross);
		fprintf( fp_dat, "%24.20e\n", h_cross);

	} 

	fclose( fp );
	fclose( fp_dat );

	gsl_odeiv_evolve_free (solver_evolve);
	gsl_odeiv_control_free (solver_control);
	gsl_odeiv_step_free (solver_step);
	return 0;
}

