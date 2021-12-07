/* $Id: conservative_dynamics.c,v 1.8 2009/03/21 18:18:47 pzimmerman Exp $ */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>

//#include "conservative_dynamics.h"

double cosu_factor( double e, double u );
double beta_phi( double e_phi );

double cosu_factor(double e, double u) 
{
	return (e * cos(u) - 1.);
}

/* Hinder et al. Eqn (A20) */
double beta_phi( double e_phi )
{
	return ( (1.0 - sqrt(1.0 - e_phi*e_phi)) / e_phi );
}
	

/* 
 * 0pN corrections to the conservative dynamics 
 */


/* Hinder et al. Eqn (A11) */
double rel_dphi_dt_0pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);
	return ( sqrt(1.-e*e) / ((u_factor*u_factor)) );
}                  	

/* Hinder et al. Eqn (A6) */
double rel_sep_0pn( double e, double u )
{
	return ( 1.0 - e*cos(u) );
}


/* 
 * 1pN corrections to the conservative dynamics 
 */


/* Hinder et al. Eqn (A22) */
double angular_ecc_1pn( double e, double eta )
{
	return ( -e * (eta - 4.) );
}

/* Hinder et al. Eqn (A2) */
double mean_motion_1pn( double e )
{
	return ( 3. / (e*e - 1.) );
}

/* Hinder et al. Eqn (A12) */
double rel_dphi_dt_1pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);

  return ( - (e*(eta - 4.) * u_factor) /
		(sqrt(1. - e*e) * u_factor*u_factor*u_factor) );
}

/* Hinder et al. Eqn (A7) */
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


/* Hinder et al. Eqn (A3) */
double mean_motion_2pn( double e, double eta )
{
	return ( ((26.*eta - 51.)*e*e + 28. - 18.)
			/ (4.*(e*e-1.) * (e*e-1.)) );
}

/* Hinder et al. Eqn (A23) */
double angular_ecc_2pn( double e, double eta )
{
	return ( (e/(96.*(e*e-1.))) * ( (41.*eta*eta-659.*eta+1152.)*e*e +
		4.*eta*eta + 68.*eta + sqrt(1.-e*e)*(288.*eta-720.) - 1248.) );
}

/* Hinder et al. Eqn (A17) */
double mean_anomaly_2pn( double e, double u, double eta, double e_phi )
{
	double u_factor = cosu_factor(e, u);
	double e_factor = 1.-e*e;
	double v_minus_u;

	v_minus_u = 2.0 * atan( (sin(u) * beta_phi(e_phi)) /
			(1.0 - beta_phi(e_phi)*cos(u)) );
	//double u_minus_v = -v_minus_u;

	return ( (1./ (8.* sqrt(e_factor) * u_factor)) * 
					( -12.* (2.* eta - 5.) * (-v_minus_u) * u_factor - 
					 e * sqrt(e_factor) * (eta - 15.) * eta * sin(u) ));
}

/* Hinder et al. Eqn (A13) */
double rel_dphi_dt_2pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);

	return ( (1. / (12.* pow(1. - e*e, 3./2.) * pow(u_factor, 5.))) *
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

/* Hinder et al. Eqn (A8) */
double rel_sep_2pn(double e, double u, double eta)
{
	double e_factor = 1.0-e*e;

	return ( (1.0 / (e_factor*e_factor)) * (1.0/72.0) * (
				 ( (8.* eta*eta + 30.*eta + 72.)*e*e + (-16.*eta*eta - 876.*eta + 756.) ) * (e*e) 
				 + (8.*eta*eta+198.*eta+360.) +
				 ( ( (-35.*eta*eta + 231.*eta - 72.)*e*e + (70.*eta*eta - 150.*eta - 468.) )*(e*e*e) + 
				 	 (-35.*eta*eta+567.*eta-648.)*e )*cos(u) + 
 			 sqrt(e_factor) * ( (360.-144.*eta*eta)*e*e + (144.*eta-360.) + 
 			 	 ( (180. - 72.*eta)*e*e + (72.*eta - 180.) )*(e*cos(u)) ) ) );
}


/* 
 * 3pN corrections to the conservative dynamics 
 */


/* Hinder et al. Eqn (A4) */
double mean_motion_3pn( double e, double eta )
{

	return ( (-1./(128.*pow(1.-e*e,7./2.))) * 
			((1536.*eta*eta*eta-3840.) *e*e*e*e + (1920.-786.*eta) *e*e
			 -768.*eta + sqrt(1.-e*e) * 
			 ( (1040.*eta*eta-1760.*eta+2496.)*e*e*e*e
			 	 + (5120.*eta*eta+123.*eta*M_PI*M_PI-17856.*eta+8544)*e*e
			 	 +896.*eta-14624.*eta+492.*eta*M_PI*M_PI-192.) + 1920.));
}

/* Hinder et al. Eqn (A24) */
double angular_ecc_3pn( double e, double eta )
{
	return ( - (e / (26880.* pow(1. - e*e, 5./2.))) *
		( ( (13440.* eta*eta + 48340.*eta - 940800.)*e*e +
		 (255360.* eta*eta + 17220.* M_PI*M_PI*eta - 28810640.* eta
		 	+ 2688000) ) * (e*e) + (-268800.*eta + 2396800.)*eta 
		  + sqrt(1. - e*e) * ( ( (1050.* eta*eta*eta -
		  		134050.*eta*eta + 786310.* eta - 8601060.)*e*e
		  	+ (-18900.* eta*eta*eta + 553980.* eta*eta + 4305.*M_PI*M_PI*eta
		  		-1246368.*eta + 2042880.) )*(e*e) + (276640.*eta +
		  	2674480. -17220.* M_PI*M_PI)*eta - 1451520.) - 17220.*M_PI*M_PI*eta - 1747200.) );
}

/* Hinder et al. Eqn (A18) */
double mean_anomaly_3pn( double e, double u, double eta, double e_phi )
{
	double u_factor = cosu_factor(e, u);
	double e_factor = 1.0-e*e;
	//double pi_factor = M_PI*M_PI*eta;
	double v_minus_u = 2.0 * atan( (sin(u) * beta_phi(e_phi)) /
			(1.0 - beta_phi(e_phi) * cos(u)) );
	double u_minus_v = -v_minus_u;

	return ( (1. / (6720. * e_factor * sqrt(e_factor) * u_factor*u_factor*u_factor)) * 
			(  35.* ( 96.* (11.*eta *eta-29.*eta+30.) * e*e 
			+ 960.*eta*eta + eta*(-13184. + 123.* M_PI*M_PI) 
			+ 8640. ) * u_minus_v * (u_factor*u_factor*u_factor) + 
			3360.* ( -12.*(2.* eta - 5.) * (u_minus_v) + 12.*e * (2.* eta - 5.) * cos(u)*u_minus_v
				+ e * sqrt(e_factor) * (eta - 15.) * eta * sin(u) ) * u_factor*u_factor +
			e * sqrt(e_factor) * (140.* ((13.*e*e-11.)*e*e - 2.)*eta*eta*eta - 
				140.*((73.*e*e -325.)*e*e + 444.)*eta*eta + 
				((3220.*e*e - 148960.)*e*e - 4305.*M_PI*M_PI + 143868.)*eta + 
				e*e * (1820.*(e*e-1.)*eta*eta*eta - 140.*(83.*e*e + 109.)*eta*eta 
					- (1120.*e*e + 4305.*M_PI*M_PI + 752.)*eta + 67200.)*cos(u)*cos(u) 
				- 2.*e*(1960.*(e*e-1.)*eta*eta*eta+6720.*(e*e-5.)*eta*eta 
					+ (-71820.*e*e - 4305.*M_PI*M_PI + 69948.)*eta + 67200.)*cos(u)
						+ 67200.)*sin(u) ) );
}


/* Hinder et al. Eqn (A14) */
double rel_dphi_dt_3pn( double e, double u, double eta )
{
	double u_factor = cosu_factor(e, u);

	return ( (1./ ( 13440.*pow(1.-e*e, 3./2.) * pow( u_factor, 7.) )) *
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
       	12915.*M_PI*M_PI*eta + 692928.*eta - 194880.) * pow(e, 7.) + 
       (4480.*eta*eta*eta + 8960.*eta*eta - 43050.*M_PI*M_PI*eta +
       	1127280.*eta - 147840.) * e*e*e*e*e) * pow( cos(u), 5.) +
      ((-16240.*eta*eta*eta - 12880.*eta*eta 
      	- 18480.*eta) * pow(e, 10.) +
       (16240.*eta*eta*eta - 91840.*eta*eta + 
       	17220.*M_PI*M_PI*eta - 652192.*eta + 100800.) * pow(e, 8.) +
       (-55440.*eta*eta*eta + 34160.*eta*eta - 
       	30135.*M_PI*M_PI*eta - 2185040.*eta + 2493120.) * e*e*e*e*e*e +
       (21480.*eta*eta*eta + 86800.*eta*eta + 163590.*M_PI*M_PI*eta
       	- 5713888.*eta + 228480.) * e*e*e*e ) * pow(cos(u), 4.) + 
      ((560.*eta*eta*eta - 137200.*eta*eta + 
      	388640.*eta + 241920.) * pow(e, 9.) +
       (30800.*eta*eta*eta - 264880.*eta*eta - 
       	68880.*M_PI*M_PI*eta + 624128.*eta + 766080.) * pow(e, 7.) +
       (66640.*eta*eta*eta + 612080.*eta*eta - 8610.*M_PI*M_PI*eta + 
       	6666080.*eta - 6652800.) * e*e*e*e*e + 
       (-30800.*eta*eta*eta-364880.*eta*eta-
       	223860.*M_PI*M_PI*eta + 9386432.*eta) * e*e*e) * pow(cos(u), 3.) +
      67200.*eta*eta + ((4480.*eta*eta*eta - 
      		20160.*eta*eta+16800.*eta) * pow(e,10.) +
      		(3920.*eta*eta*eta+475440.*eta*eta - 
      		 17220.*M_PI*M_PI*eta + + 831952.*eta - 7257600.) * pow(e,8.) + 
      		(-75600.*eta*eta*eta + 96880.*eta*eta + 
      		 154980.*M_PI*M_PI*eta - 
      		 3249488.*eta-685440.) * e*e*e*e*e*e +
      		(5040.*eta*eta*eta - 659120.*eta*eta + 
      		 25830.*M_PI*M_PI*eta - 7356624.*eta + 6948480.) * e*e*e*e
      		+ (-5040.*eta*eta*eta + 190960.*eta*eta + 
      			137760.*M_PI*M_PI*eta - 
      			7307920.*eta + 107520.) * e*e ) * pow(cos(u), 2.) -
      761600.*eta + ((-2240.*eta*eta*eta - 168000.*eta*eta - 
      	424480.*eta) * pow(e, 9.) + (28560.*eta*eta*eta + 
      		242480.*eta*eta + 34440.* M_PI*M_PI*eta - 1340224.* eta +
      		725760.) * pow(e, 7.) + (-33040.* eta*eta*eta-754880.* eta*eta
      			-172200.*M_PI*M_PI*eta+5458480.*eta-221760.)* e*e*e*e*e + 
      		(40880.*eta*eta*eta+738640.*eta*eta + 30135.*M_PI*M_PI*eta
      		 + 1554048.*eta - 2936640.) *e*e*e + (-560.*eta*eta*eta -
      		 	 100240.*eta*eta - 43050.*M_PI*M_PI*eta + 
      		 	 3284816.*eta - 389760.) * e ) * cos(u) + 
      sqrt(1.-e*e) * ( ((-127680.*eta*eta + 544320.*eta-739200.) * pow(e,7.) +
      			(-53760.*eta*eta - 8610.*M_PI*M_PI*eta + 674240.*eta - 
      		 	 67200.) * e*e*e*e*e ) * pow(cos(u), 5.) + ( (161280.*eta*eta - 
      		 	 	 477120.*eta + 537600.) * pow(e,8.) + (477120.*eta*eta + 
      		 	 	 	17220.*M_PI*M_PI*eta - 2894080.*eta + 2217600.) * pow(e,6.)
      		 	 + (268800.* eta*eta + 25830.* M_PI*M_PI*eta - 2721600.*eta +
      		 	 	 1276800.) * e*e*e*e ) * pow(cos(u),4.) + 
      		((-524160.*eta*eta + 1122240.*eta - 940800.) * pow(e,7.) 
      		 + (-873600.*eta*eta - 68880.* M_PI*M_PI*eta + 7705600.*eta
      		 	 - 3897600.) * e*e*e*e*e + (-416640.*eta*eta - 
      		 	 	 17220.* M_PI*M_PI*eta + 3357760.*eta - 
      		 	 	 3225600.) * e*e*e ) * cos(u)*cos(u)*cos(u) + 
      		((604800.*eta*eta - 50400.*eta - 403200.) * pow(e,6.) 
      		 + (1034880.*eta*eta + 103320.*M_PI*M_PI*eta 
      		 	 - 11195520.*eta + 5779200.) * e*e*e*e +
      		 (174720.*eta*eta - 17229.*M_PI*M_PI
      		 	- 486080.*eta + 2688000.) * e*e ) * cos(u)*cos(u)
      		+ ((-282240.*eta*eta - 450240.*eta + 1478400.) *e*e*e*e*e
      			+ (-719040.*eta*eta - 68880*M_PI*M_PI*eta
      				+ 8128960.*eta - 5040000.) * e*e*e +
      			(94080.*eta*eta + 25830.*M_PI*M_PI*eta 
      			 - 1585920.*eta - 470400.) * e ) * cos(u) 
      		- 67200.*eta*eta + 761600.*eta + 
      		e*e*e*e * (40320.*eta*eta + 309120.*eta -672000.) + e*e * (208320.*eta*eta + 
      					17220.* M_PI*M_PI*eta - 2289280.*eta + 
      					1680000.) - 8610.*M_PI*M_PI*eta - 201600.)
      		+ 8610.* M_PI*M_PI*eta + 201600.) );
}

/* Hinder et al. Eqn (A9) */
double rel_sep_3pn( double e, double u, double eta )
{
	double pi_factor = M_PI*M_PI*eta;
	double e_factor = 1.0 - e*e;

	return ( (1. / (181440.*sqrt(e_factor)*e_factor*e_factor)) *
	 ( ( (-665280.*eta*eta + 1753920.*eta - 1814400.) * e*e*e*e +
			 (725760.*eta*eta - 77490.*pi_factor + 5523840.*eta - 3628800.) * e*e +
			 (544320.*eta*eta + 154980.*pi_factor - 14132160.*eta + 7257600.) ) * e*e 
		- 604800.*eta*eta + 6854400.*eta + 
		 ( ( (302400.*eta*eta - 1254960.*eta + 453600.) * e*e*e*e +
		 	 (-1542240.*eta*eta - 38745.*pi_factor + 6980400.*eta - 453600.) * e*e +
		 	 (2177280.*eta*eta + 77490.*pi_factor - 12373200.*eta + 4989600.) ) * e*e*e + 
		 	(-937440.*eta*eta*eta-37845.* pi_factor + 6647760.*eta-4989600.) * e ) * cos(u) +
		 sqrt(e_factor) * ( ( (-4480.*eta*eta*eta-25200.*eta*eta+22680.*eta-120960.) * e*e*e*e +
		 		 (13440.*eta*eta*eta+4404960.-eta*eta+116235*pi_factor-12718296.*eta+5261760.) * e*e +
		 		 (-13440.*eta*eta*eta + 2242800.*eta*eta+348705.*pi_factor-19225080.*eta+1614160.) ) * e*e
		 	 + 4480.*eta*eta*eta+45360.*eta*eta-8600904.*eta + 
		 	 ( ( (-6860.*eta*eta*eta - 550620.*eta*eta-986580.*eta + 120960.) * e*e*e*e +
		 	 		(20580.*eta*eta*eta - 2458260.*eta*eta + 3458700.*eta - 2358720.) * e*e +
		 	 		(-20580.*eta*eta*eta-3539340.*eta*eta-116235.*pi_factor+20173860.*eta-16148160.) ) * e*e*e
		 	 	 + (6860.*eta*eta*eta-1220940.*eta*eta+464940.*pi_factor+17875620.*eta-417440.) * e ) * cos(u)
       + 116235.*pi_factor + 1814400.) - 77490.*pi_factor - 1814400. ) );
}
