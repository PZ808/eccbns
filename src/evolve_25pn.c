/*
 * $Id: evolve_25pn.c,v 1.40 2009/02/21 20:06:38 pzimmerman Exp $ 
 *
 * evolve_25pn.c
 *
 * Program numerically evolves a system of ODE's 
 * and generates a 2.5-pN eccentric binary waveform. 
 * We refer to arxiv:gr-qc/0712.3199v1 
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include "mean_motion_init.h"
#include "mass_params.h"

/*
 * function: eccentric_binary_evolution()
 * 
 * Evolution Equations for the mean-motion and eccentricity, n(t) and e(t).
 * 		- 2pN (beyond RR order) Adiabatic Evolution Equations in the Harmonic Gauge
 *		- Secular contributions are considered
 *  	- Oscillatory contibutions are ignored 
 * 		- Tail effects arising at the 1.5 pN order are not treated
 *	 	- Geometrized units are used
 */


int eccentric_binary_evolution (double t, const double y[], double f[],
		void* params)
{
	struct mass_params *mp = (struct mass_params *) params;  
	double M   = mp->M;
	double eta = mp->eta;
	int PNorder = mp->PNorder; 

	double n   = y[0];   
	double e   = y[1];   

	double dn_dt_0pn, dn_dt_1pn, dn_dt_2pn;
	double de_dt_0pn, de_dt_1pn, de_dt_2pn;

	double xi = M * n;
	double e_factor = 1. - e*e;
	double dn_dt_factor = pow(xi, 5./3.) * n*n * eta;
	double de_dt_factor = -pow(xi, 5./3.) * n * eta * e;

	dn_dt_0pn =  (96. + 292.* e*e + 37.* pow(e,4))  / ( 5.* pow(e_factor, 7./2.));

	dn_dt_1pn = ( pow(xi, 2./3.) / (280.* pow(e_factor,9./2.) )) * (20368. - 14784.*eta + (219880. - 159600.*eta) *e*e 
			+ (197022. - 141708.*eta) * pow(e,4) + (11717. - 8288.*eta) * pow(e,6) );  

	dn_dt_2pn = ( pow(xi, 4./3.) / (30240. * pow(e_factor, 11./2.)) ) * ( 12592864. - 13677408.*eta + 1903104.*eta*eta
			+ (133049696. - 185538528.*eta + 61282032.*eta*eta) * e*e + (284496744. - 411892776.*eta + 166506060.*eta*eta) * pow(e,4)
			+ (112598442. - 142089066.*eta + 64828848.*eta*eta) * pow(e,6) + (3523113. - 3259980.*eta + 1964256.*eta*eta) * pow(e,8)
			+ 3024.* (96. + 4268.* e*e + 4386.* pow(e,4) + 175.*pow(e,6)) * (5. - 2.*eta) * sqrt(e_factor) );


	de_dt_0pn = (304. + 121.* e*e) / (15.* pow( e_factor, 5./2.));

  de_dt_1pn = ( pow(xi, 2./3.) / (2520.* pow(e_factor,7./2.)) ) * ( 340968. - 228704.*eta + (880632. - 651252.*eta) *e*e 
			+ (125361. - 93184.*eta) * pow(e,4) );

	de_dt_2pn = ( pow(xi, 4./3.) / (30240.* pow(e_factor, 9./2.)) ) * (20815216. - 25375248.*eta + 4548096.*eta*eta 
			+ (87568332. - 128909916.*eta + 48711348.*eta*eta) *e*e + (69916862. - 93522570.*eta + 42810096.*eta*eta) * pow(e,4)
			+ (3786543. - 4344852.*eta + 2758560.*eta*eta) * pow(e,6) + 1008.* (2672. + 6963.*e*e + 565.*pow(e,4)) * (5. - 2.*eta) 
			* sqrt(e_factor) );


  if ( PNorder == 0 ) 
	{
		f[0] = dn_dt_factor * dn_dt_0pn; 
		f[1] = de_dt_factor * de_dt_0pn;
	} 
	else if ( PNorder == 1 ) 
	{  
    f[0] = dn_dt_factor * (dn_dt_0pn + dn_dt_1pn);
    f[1] = de_dt_factor * (de_dt_0pn + de_dt_1pn);
	}
	else if ( PNorder == 2 )  
	{
    f[0] = dn_dt_factor * (dn_dt_0pn + dn_dt_1pn + dn_dt_2pn);
    f[1] = de_dt_factor * (de_dt_0pn + de_dt_1pn + de_dt_2pn);
	}
	else 
	{
		fprintf( stderr, "Error in post-Newtonian order: %d\n", PNorder );
		exit(1);
	}  

	// fprintf( stderr, "dn_dt = %.6f\n", f[0] );
	// fprintf( stderr, "de_dt = %.6f\n", f[1] );

	return GSL_SUCCESS;
}

int main ( int argc, char *argv[] )
{


	/*
	 * input and output variables
	 */

	/* input binary parameters */
	double M;     			/* total mass            */
	double eta;   			/* symmetric mass ratio  */
	double f_init;      /* init. gw freq. 		   */
	double e;     			/* eccentricity        	 */
	double p;     			/* semi-latus rectum   	 */

	int pn_order; 

	double xi;    /* dimensionless parameter determined by M and n       */
	double chi;   /* dimensionless factor which simplifies the equations */

	/* output file parameters */
	FILE* fp = NULL;
	FILE* fp_dat = NULL;
	const int fname_len = 256;
	char fname[fname_len];

	/*
	 *	variables needed to integrate equations of motion
	 */

	/* time stepping variables for integrating ode's */
	int i;
	double srate = 16384.0;        /* sample rate (Hz)                    */
	double dt = 1.0/srate;         /* sampling interval (s)               */
	double t, t_next;              /* time index variables                */
	double step, step_size;        /* integrator step size and result     */
	double n_final;                /* final mean motion                   */

	/* dynamical equations for integration */
	double y[2];    

	/* conversion to and from dimensionless units */
	const double c_si = 299792458;        /* m s^{-1}                     */
	const double g_si = 6.67259e-11;      /* kg^{-1} m^3 s^{-2}           */
	const double t_sun = 4.92549095e-6;   /* mass of the sun in seconds   */
	const double m_sun = 1.98892e30;      /* kg                           */
  /* distance normalization */
	const double R = 1.0e-30;  			 
	const double iota = 0.25 * M_PI;

	/*
	 * variables needed to solve kepler's equation 
	 */

	/* variables needed for Mikkola's method */
	double a, b, sgn_b, z, s;
	/* mean anomaly: l = f(n) + pn corrections */
	double l, sgn_l;
	double l_2pn;
	/* dimmensionless parameter introduced to simplify things */
	double beta;
	/* true anomaly: v - u = 2 * atan( (beta * sin( u ) ) / ( 1 - beta * cos(u) ) )  */
	double v_minus_u;
	/* eccentric anomaly */
	double u;

	/*
	 * 	dynamical variables required to generate GW polarizations 
	 */

	double r, phi;
	double phi_0pn;
  double r_0pn, r_1pn, r_2pn;                		
	double dr_dt, dphi_dt;
	double dr_dt_0pn, dr_dt_1pn, dr_dt_2pn,
				 dphi_dt_0pn, dphi_dt_1pn, dphi_dt_2pn;
	double lambda;
	double k;
	double W;
	double W_0pn, W_1pn, W_2pn;
	double h_plus;          /* plus polarizaion  */
	double h_cross;         /* cross polarizaion */
	//double gw_cycles; 

	struct mass_params input_params;

	/*
	 * GSL differential equation solving mechanism 
	 * Solver step-type: GSL Embedded Runge-Kutta-Fehlberg (4,5) 
	 * Solver control: GSL standard adaptive control
	 */

	const gsl_odeiv_step_type* solver_type 
		= gsl_odeiv_step_rkf45;
	gsl_odeiv_step* solver_step 
		= gsl_odeiv_step_alloc( solver_type, 2 );
	gsl_odeiv_control* solver_control
		= gsl_odeiv_control_standard_new( 0.0, 1.0e-10, 1.0, 1.0 );
	gsl_odeiv_evolve* solver_evolve
		= gsl_odeiv_evolve_alloc( 2 );
	gsl_odeiv_system solver_system = { eccentric_binary_evolution,
		NULL, 2, &input_params };

	/* 
	 * parse the command line arguments and open the output file
	 */

	if ( argc != 6 )
	{
		fprintf( stderr, "error: incorrect number of arguments\n" );
		fprintf( stderr, "usage: %s M eta f_init e pn_order\n", argv[0] );
		return 1;
	} 

	M        = atof( argv[1] ); 	/* mass [solar masses]        			*/
	eta      = atof( argv[2] ); 	/* symmetric mass ratio [dimensionless] */
	f_init   = atof( argv[3] );   /* starting frequency [Hz]			 	*/
	e        = atof( argv[4] );   /* eccentricity [dimensionless] 		*/
	pn_order = atoi( argv[5] );

  input_params.M = M;
  input_params.eta = eta;
  input_params.PNorder = pn_order;

	fprintf( stderr, "M = %4.2f solar mass \n", M );
	fprintf( stderr, "eta = %.6f\n", eta );
	fprintf( stderr, "f_init = %4.2f Hz \n", f_init );
	fprintf( stderr, "e = %4.2f\n", e );
	fprintf( stderr, "pN_order above RR = %2d\n", pn_order );

	/* create the file name and file pointer */
	snprintf( fname, fname_len * sizeof(char), 
			"evolve_%gpN_%4.2f_%4.2f_%4.2f_%4.2f.txt", 
			(float)pn_order, M, f_init, eta, e );
	fp = fopen( fname, "w" );

#if 0
	snprintf( fname, fname_len * sizeof(char), 
			"evolve_%2.1fpN_%4.2f_%4.2f_%4.2f_%4.2f.dat", 
			(float)pn_order, M, f_init, eta, e );
	fp_dat = fopen( fname, "w" );
	fprintf( fp_dat, "%% dx = %32.16e\n", dt );
#endif

	/*
	 * Make the root solving mechanism
	 */

	int n_init_status;
	int iter = 0;
	int max_iter = 1000;
	const gsl_root_fsolver_type *root_finder_type;
	gsl_root_fsolver *root_solver;
	double n = 0.;
  double n_lo = 100.0;       	/* estimated lower bracket */
  double n_hi  = 200.0;     	/* estimated upper bracket */
	double root_eps_abs = 0.0;
	double root_eps_rel = 1.0e-6;

  gsl_function F;
	struct ic_params n_params = {f_init, M, eta, e, pn_order};

  F.function = &n_init;
	F.params   = &n_params;

	root_finder_type = gsl_root_fsolver_brent;
	root_solver = gsl_root_fsolver_alloc (root_finder_type);  
	gsl_root_fsolver_set( root_solver, &F, n_lo, n_hi );

	iter = 0;
  do {
    n_init_status = gsl_root_fsolver_iterate (root_solver);
		n = gsl_root_fsolver_root (root_solver);
		n_lo = gsl_root_fsolver_x_lower (root_solver);
		n_hi = gsl_root_fsolver_x_upper (root_solver);
    n_init_status = gsl_root_test_interval ( n_lo, n_hi,
				root_eps_abs, root_eps_rel);

		if ( n_init_status == GSL_SUCCESS )
		{
			fprintf( stderr, "n_initial = %e \n", n*t_sun );
			break;
		}

		if ( iter == max_iter || n_init_status != GSL_CONTINUE ) 
		{
			fprintf( stderr, " Error in GSL root finder \n" );
			exit(1);
		}
	} while (1);			

	/*
	 * 	set initial values
	 */

	t = 0.0;

	y[0] = n * t_sun;
	y[1] = e;

  /* termination condition: Schwarzschild Radius */
	n_final = 1. / ( pow(6.0, 3./2.) * M );

	/*
	 * generate waveform
	 */

	i = 0;
	step = dt / t_sun;


	while ( 1 )
	{
		/* advance the time */ 
		t = i * step;   
		t_next = ++i * step;
		step_size = step;

		//fprintf( stderr, "y[0] is %f, y[1] is %f\n", y[0], y[1] );

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
			//fprintf( stderr, "N_gw_cycles = %e\n", phi / M_PI);
			fprintf( stderr, "breaking due to (n = %e) > (n_final = %e)\n", y[0], n_final );
			break;       
		}

    p = ( 1. - y[1]*y[1] ) / pow( y[0] * M, 2./3.);


		/*
		 * Use Minkkola's method to solve Kepler's equation
		 */

		/* compute the mean anomaly from the mean motion */
		//l = y[0] * t;

    //fprintf( stderr, "mean_anomaly is %e\n", l);

		/* 2pN accurate Kepler Equation in the harmonic gauge: equation (13) G.T */
		l_2pn = u - y[0] * sin(u) +
			( pow(xi, 4./3.) / (8.* sqrt(1. - y[1]*y[1]) * chi) ) 
			* ( (15.* eta - eta*eta) * y[1]*sin(u) * sqrt(1. - y[1]*y[1])
					+ 12.* (5. - 2.* eta) * v_minus_u * chi ); 

		/* range reduction of l */

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
			sgn_l = 1.0;
		else
		{
			sgn_l = -1.0;
		}

		l *= sgn_l;

		/* compute alpha and beta of Minkkola Eq. (9a) */
		a  = ( 1.- y[1] ) / ( 4.* y[1] + 0.5);
		b  = ( 0.5 * l ) / ( 4.* y[1] + 0.5 );

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
		s = s - 0.078 * pow( s, 5. ) / ( 1.+ y[1] );

		/* Finally Mikkola Eq. (8) gives u */
		u = l + y[1] * ( 3.* s - 4.* pow( s, 3. ) ) ;   

    //fprintf( stderr, "u = %e\n", u);

		/* correct the sign of u */
		u *= sgn_l;

		/*
		 * compute the variables needed to get the waveform
		 */

		/* beta blows up for zero eccentricity */

		beta = ( 1. - sqrt(1. - y[1]*y[1]) ) / y[1];
		v_minus_u =  2.* atan( ( beta * sin(u) ) / 
				( 1. - beta * cos(u) ) );  

		xi = M * y[0];  
		chi = ( 1. - y[1] * cos(u) );

		k = ( 3.* pow(xi, 2./3.) ) / ( 1. - y[1]*y[1] ) +
			( pow(xi, 4./3.) / ( 4.* pow( 1.* - y[1]*y[1], 2.) ) ) * 
			( 78. - 28.* eta + (51. - 26.* eta) * y[1]*y[1] );

    lambda = ( 1.0 + k ) * l;

		//fprintf( stderr, "beta = %e, v_minus_u = %e, k = %e lamda = %e\n", beta, v_minus_u, k, lambda );	
    if ( isinf( beta ) || isinf( v_minus_u ) || isnan( k ) || isnan( lambda ) )
    {
      fprintf( stderr, "beta is %f, v_minus is %f, k is %f, lambda is %f\n", beta, v_minus_u, k, lambda );
      exit ( 1 );
    }

		/*
		 * 	variables at the 0pN level 
		 */

		W_0pn  = v_minus_u + y[1] * sin(u);

		r_0pn  		= pow( ( g_si * M * m_sun ) / ( y[0]*y[0] ), 1./3. ) * chi;
		dr_dt_0pn 	= ( pow( g_si * M * m_sun * y[0], 1./3. ) * (y[1] * sin(u)) ) / chi;

		phi_0pn		= l + W_0pn;

		dphi_dt_0pn = ( y[0] * sqrt(1. - y[1]*y[1]) ) / ( chi*chi ); 

		//fprintf( stderr, " W_0pn = %e, phi_0pn = %e\n",  W_0pn, phi_0pn );	

		/*
		 * 	variables at the 1pN level
		 */

		W_1pn = ( ( 3.* pow(xi, 2./3.) ) / ( 1. - y[1]*y[1] ) ) * ( v_minus_u + y[1] * sin(u) );

		r_1pn = r_0pn * ( pow(xi, 2./3.) 
		   	* (-18. + 2.* eta - (6. - 7.* eta) * y[1]*cos(u)) ) / (6.* chi); 

		dr_dt_1pn 	= dr_dt_0pn * ( (pow(xi, 2./3.) / 6.) * (6. - 7.* eta) );

    dphi_dt_1pn = dphi_dt_0pn * ( pow(xi, 2./3.) / ( 32.* (1. - y[1]*y[1]) * chi) ) *
			( 3. - (4. - eta) * y[1]*y[1] + (1. - eta) * y[1]*cos(u) );

		/*
		 * 	variables at the 2pN level
		 */

		W_2pn = pow(xi, 4./3.) / ( 32.* pow(1. - y[1]*y[1], 2.) * pow(chi, 3.) ) 
			* ( 8.* ( 78. - 28.* eta + (51. - 26.* eta) * y[1]*y[1] 
						- 6.* (5. - 2.*eta) * pow(1. - y[1]*y[1], 3./2.) ) * v_minus_u * pow(chi, 3.) );

		r_2pn = r_0pn * ( pow(xi, 4./3.) / (72.* (1. - y[1]*y[1]) * chi) ) * 
		  ( -72.* (4. - 7.* eta) + ( 72. + eta * (30. + 8.* eta)
		   		 											 - (72. - eta * (231. + 35.* eta) * (y[1] * cos(u))) ) 
				- 36.* (5. - 2.* eta) * (2. + y[1] * cos(u)) * sqrt(1. - y[1]*y[1]) );

		dr_dt_2pn =  dr_dt_0pn * ( pow(xi, 4./3.) / ( 72.* pow(chi, 3. )) ) *
			( -468. - 15.*eta + 35.* eta*eta + eta * (135. - 9.*eta) * y[1]*y[1] 
				+ ( 324. + eta * (342. - 96.* eta) * y[1] * cos(u) 
					+ (216. - 693.* eta + 105.* eta*eta) * pow(y[1]*cos(u), 2.) 
					- (72. - 231.* eta + 35.* eta*eta) * pow(y[1]*cos(u), 3.) 
					+ ((36.* chi*chi * (4. - y[1]*cos(u)) * (5. - 2* eta)) / sqrt(1. - y[1]*y[1])) ) );

		dphi_dt_2pn = dphi_dt_0pn * 
			( pow(xi, 4./3.) / ( 12.* ( 1. - y[1]*y[1]) * pow(chi, 3.) ) ) *
			( 144. - 48.* eta - (162. + 68.* eta -2.* eta*eta) * y[1]*y[1] 
				+ (60. + 26.* eta - 20.* eta*eta) * pow( y[1], 4.) + (18.* eta +12.* eta*eta) * pow(eta, 6.) 
				+ ( -216. + 125.*eta + eta*eta + (102. + 188.* eta + 16.* eta*eta) * y[1]*y[1] 
					- (12.* + 97.* eta - eta*eta) * pow(y[1], 4.) ) * y[1]*cos(u) 
  			+ (108. - 97.* eta - 5.*eta*eta + (66. - 136.*eta + 4.* eta*eta) * y[1]*y[1] 
					- (48. - 17.* eta + 17.* eta*eta) * pow(y[1], 4.)) * pow( y[1]*cos(u), 2.) 
				+ ( -36. + 2.* eta - 8.* eta*eta - (6. - 70.* eta - 14.* eta*eta) * y[1]*y[1] ) * pow(y[1]*cos(u), 3.) 
				+ 18.* pow(chi, 2.) * ( 1 - 2.* y[1]*y[1] + y[1]*cos(u)) * (5. - 2.* eta) * sqrt(1. - y[1]*y[1]) ); 

		if ( pn_order == 0 )
		{
	    r = r_0pn;
			phi	= phi_0pn;
			dr_dt   = dr_dt_0pn;
			dphi_dt = dphi_dt_0pn;
    }
		else if ( pn_order == 1 )
		{
			W = W_0pn + W_1pn;
			r = r_0pn + r_1pn;
			dr_dt   = dr_dt_0pn + dr_dt_1pn;
			phi	    = lambda + W; 
			dphi_dt = dphi_dt_0pn + dphi_dt_1pn;
		}	
		else if ( pn_order == 2 )
		{

			W 		= W_0pn + W_1pn + W_2pn;
			r	   	= r_0pn + r_1pn + r_2pn;
			dr_dt   = dr_dt_0pn + dr_dt_1pn + dr_dt_2pn;
			phi	    = lambda + W;
			dphi_dt = dphi_dt_0pn + dphi_dt_1pn + dphi_dt_2pn; 
    }
		else
		{
			fprintf( stderr, "error in pN order" );
			exit( 1 );
		} 		

		h_plus 	= - ( eta * M / R ) * ( (1.+ cos(iota)*cos(iota)) *
				((M/r + r*r * dphi_dt*dphi_dt - dr_dt*dr_dt) * cos(2.*phi) + 
				 2.* r * dr_dt*dphi_dt * sin(2.* phi)) +  
				sin(iota)*sin(iota) * (M/r - r*r * dphi_dt*dphi_dt - dr_dt*dr_dt) );                  

#if 0
		h_cross	= - ( (2.* g_si * M * eta) / (c_si*c_si*c_si*c_si * R) )
			* ( ((g_si * M / r) + r*r * dphi_dt*dphi_dt - dr_dt*dr_dt) * sin(2.* phi) - 
			   	2.* dr_dt * r * dphi_dt * cos(2.* phi) );

#endif
#if 0
		fprintf( stderr, "t=%f, f_gw=%e, n=%f, e=%f, a=%e, b=%e, u=%e\n",
				t, y[0] / M_PI, y[0], y[1], a, b, u );
		fprintf( fp, "%e %e %e %e %e\n", t, y[0], y[1], h_plus, h_cross);
#endif
		fprintf( fp, "%32.16e %32.16e %32.16e %36.16e\n", t_sun * t, y[0], p, h_plus);

		/* fprintf( fp, "%32.16e %32.16e %32.16e %36.16e %32.16e\n", t, y[0], y[1], p, h_cross); */
		/* fprintf( fp_dat, "%32.16e\n", h_cross); */
	} 

	fclose( fp );
#if 0
	fclose( fp_dat );
#endif
	gsl_root_fsolver_free (root_solver); 
	gsl_odeiv_evolve_free (solver_evolve);
	gsl_odeiv_control_free (solver_control);
	gsl_odeiv_step_free (solver_step);
	return 0;
}

