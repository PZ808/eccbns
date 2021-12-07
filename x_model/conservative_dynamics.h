/* $Id: conservative_dynamics.h,v 1.2 2009/04/02 13:43:28 pzimmerman Exp $ */

/* angular eccentricity *
 * EQNS (A21-A24)       */

double angular_ecc_1pn( double e, double eta );
double angular_ecc_2pn( double e, double eta );
double angular_ecc_3pn( double e, double eta );

/* mean anomaly   *
 * EQNS (A16-A18) */

double mean_anomaly_2pn( double e, double e_phi, 
		double u, double eta );
double mean_anomaly_3pn( double e, double e_phi, 
		double u, double eta );

/* mean motion  *
 * EQNS (A2-A4) */

double mean_motion_1pn( double e );
double mean_motion_2pn( double e, double eta );
double mean_motion_3pn ( double e, double eta );

/* Relative rate of * 
 * change of phase  *
 * EQNS (A11-A14)   */

double rel_dphi_dt_0pn( double e, double u, 
		double eta );
double rel_dphi_dt_1pn( double e, double u, 
		double eta );
double rel_dphi_dt_2pn( double e, double u, 
		double eta );

/* Relative rate of     * 
 * change of separation *
 * EQNS (A6-A9)         */

double rel_sep_0pn( double e, double u );
double rel_sep_1pn( double e, double u, 
		double eta );
double rel_sep_2pn( double e, double u,
		double eta );
double rel_sep_3pn( double e, double u,
		double eta );
