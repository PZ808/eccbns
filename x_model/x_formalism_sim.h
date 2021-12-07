/* $Id: x_formalism_sim.h,v 1.1 2009/04/02 13:43:28 pzimmerman Exp $ */
double mean_motion(
		int cPN,
		double eta,
		double x,
		double e );

double dphi_dt(
		int cPN,
		double u,
		double eta,
		double x,
		double e );

double separation( 
		int cPN,
		double u, 
		double eta,
		double x,
		double e );

double mikkola_finder(
		double e,
		double l );

struct ode_parameters
	{
		double eta;
		double ecc_anomaly;
		int conservative_pn_order;
		int radiation_pn_order;
	};
