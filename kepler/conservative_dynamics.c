/* 
 *	conservative_dynamics.c 
 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "conservative_dynamics.h"

double cosu_factor( double e, double u );

double cosu_factor(double e, double u) 
{
	return (e * cos(u) - 1.);
}

/* 0pN corrections to the conservative dynamics */

double rel_dphi_dt_0pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);
	return ( sqrt(1.-e*e) / ((u_factor*u_factor)) );
}                  	
double rel_sep_0pn( double e, double u )
{
	return ( 1.0 - e*cos(u) );
}

/* 1pN corrections to the conservative dynamics */

double angular_ecc_1pn( double e, double eta )
{
	return ( -e * (eta - 4.) );
}

double mean_motion_1pn( double e )
{
	return ( 3. / (e*e - 1.) );
}

double rel_dphi_dt_1pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);

  return ( - (e*(eta - 4.) * u_factor) /
		(sqrt(1. - e*e) * u_factor*u_factor*u_factor) );
}

double rel_sep_1pn( double e, double u, double eta )
{ 
	double u_factor = cosu_factor(e, u);

	return ( 2.*(u_factor) / (e*e - 1.) + 
	 	(1./6.) * (2.*(eta - 9.) +
		e * (7.*eta - 6.) * cos(u)) );
}	

/* 
 * 2pN corrections to the conservative dynamics 
 */

double mean_motion_2pn( double e, double eta )
{
	return ( ((26.*eta - 51.)*e*e + 28. - 18.)
			/ (4.*(e*e-1.) * (e*e-1.)) );
}

double angular_ecc_2pn( double e, double eta )
{
	return ( (e/(96.*(e*e-1.))) * ( (41.*eta*eta-659.*eta+1152.)*e*e +
		4.*eta*eta + 68.*eta + sqrt(1.-e*e)*(288.*eta-720.) - 1248.) );
{

double mean_anomaly_2pn( double e, double u, double eta,
		double u_minus_v(double u, double beta_phi) )
{
	double u_factor = cosu_factor(e, u);
	double e_factor = 1.-e*e;

	return ( (1./ (8.* sqrt(e_factor) * u_factor)) * 
					( -12.(2.*eta-5.)*v_minus_u*u_factor - 
					 e*sqrt(e_factor)*(eta-15.)*eta*sin(u) ));
}

double rel_dphi_dt_2pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);
	double e_factor = 1.-e*e;

	return ( (1./(12.*pow(1.-e*e, 3./2.)*pow(u_factor, 5.))) *
		( (-12.*eta*eta-18.*eta)*pow(e,6) + (20.*eta*eta-26.*eta-60.)*pow(e,4) +
		(-2.*eta*eta+50.*eta+75.)*e*e + ((-14.*eta*eta+8.*eta-147.)*pow(e,5) +
		(8.*eta*eta+22.*eta+42.)*e*e*e) * cos(u)*cos(u)*cos(u) + 
		((17.*eta*eta-17.*eta+48.)*pow(e,6) + (-4.*eta*eta-38.*eta+153.)*e*e*e*e +
		(5.*eta*eta-35.*eta+114.)*e*e) * cos(u)*cos(u) - 36.*eta + 
		((-eta*eta+97.*eta+12.)*pow(e,5) + (-16.*eta*eta-74.*eta-81.)*e*e*e +
		(-eta*eta+67.*eta-246.)*e) * cos(u) + 
		sqrt(1.-e*e) * (e*e*e*(36.*eta-90.)*cos(u)*cos(u)*cos(u) +
    ((180.-72.*eta)*e*e*e*e + (90.-36.*eta)*e*e)*cos(u)*cos(u) + 
    ((144.*eta-360.)*e*e*e + (90.-36.*eta)*e)*cos(u) +
    e*e*(180.-72.*eta) + 36.*eta - 90.) + 90.) );
}

double rel_sep_2pn(double e, double u, double eta)
{
	double u_factor = cosu_factor(e, u);

	return ( (1./(1.-e*e)*(1.-e*e)) *
			( (1./72.)*(8.*eta*eta+30.*eta+72.)*e*e*e*e +
        (1./72.)*(-16.*eta*eta-876.*eta+756.)*e*e +
 				(1./72.)*(8.*eta*eta+198.*eta+360.) + 
       ( (1./72.)*(-35.*eta*eta+231.*eta-72.)*e*e*e*e*e +
        (1./72.)*(70.*eta*eta-150.*eta-468.)*e*e*e +
 				(1./72.)*(-35.*eta*eta+567.*eta-648.)*e )*cos(u) + 
 			 sqrt(1.-e*e) * ( (1./72.)*(360.-144.*eta*eta)*e*e +
 			 	   (1./72.)*(144.*eta-360.) + 
 			  	((1./72.)*(180.-72.*eta)*e*e*e + 
 			  	 (1./72.)*(72.*eta-180.)*e)*cos(u) ) ) ):
}

/* 3pN corrections to the conservative dynamics */

double mean_anomaly_3pn( double e, double u, 
		double u_minus_v(double u, double beta_phi) )
{
	return ( (1./ (8.* sqrt(e_factor) * u_factor*u_factor*u_factor)) * 
			(35.*(96.*(11.*eta *eta-29.*eta+30.) *e*e 
			+960.*eta*eta+*eta*(-13184.+123.*M_PI*M_PI)) ) );
}


double rel_dphi_dt_3pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);

	return ( (1./ ( 13440.*pow(1.-e*e, 3./2.) * pow( u_factor, 7.) ) *
		( (10080.*eta*eta*eta + 40320.*eta*eta 
			 - 15120.*eta) * pow(e, 10.) +
			(-52640.*eta*eta*eta + 13440.*eta*eta + 
			 483280.*eta) * pow(e, 8.) +
			(84000.*eta*eta*eta - 190400.*eta*eta - 
			 17220.*M_PI*M_PI*eta - 50048.*eta-241920.) * pow(e, 6.) +
			(-52640.*eta*eta*eta + 516880.*eta*eta + 
			 68880.*M_PI*M_PI*eta - 1916048.*eta - 262080.) * e*e*e*e +
      (4480.*eta*eta*eta-412160.*eta *eta - 
       30135.*M_PI*M_PI*eta + 553008.*eta + 342720.) * e*e + 
      ((13440.*eta*eta*eta + 94640.*eta*eta - 
      	113680.*eta - 221760.) * pow(e, 9.) +
       (-11200.*eta*eta*eta - 112000.*eta*eta + 
       	12915.*M_PI*M_PI*eta + 692928.*eta - 194880.) * pow(e, 7.)
