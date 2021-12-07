/* $Id: x_ecc_ode_pn_terms.c,v 1.13 2009/05/27 05:54:11 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_bessel.h>

#include "x_drive.h"

/*
 * Abromowitz and Stegun 'Handbook of Mathematical Functions' pp 358
 * Recursion relations for Bessel Functions (9.1.27)
 *
 * J_{n+1} + J_{n-1} = \frac{2n}{x} J_n
 * J_{n-1} - J_{n+1} = 2 J_n^{'} 
 *
 * J_n^{'} = J_{n-1} - \frac{n}{x} J_n
 * J_n^{'} = -J_{n+1} + \frac{n}{x} J_n
 *
 */

static double kappa_e_part_sum( double p, double e, double J_p, double J_p_prime )
{
  double kappa_e = 0;
  double over_e = 1.0/e;
	double e_pow_2 = e*e;
  double over_ee = 1.0/e_pow_2;
  double p_pow_2 = p*p;

  kappa_e = 0.25 * p*p_pow_2 * ( 
      ( (-e_pow_2 + (-3.0 + over_ee) * over_ee + 3.0) * p_pow_2 +
        1.0/3.0 + (-1.0 + over_ee) * over_ee ) * J_p * J_p +
      ( -3.0 * e + (-4.0 * over_ee + 7.0) * over_e) * p * J_p_prime * J_p +
      ( (e_pow_2 + over_ee - 2.0) * p_pow_2 + over_ee  - 1.0 ) * J_p_prime * J_p_prime );

  return kappa_e;
}

static double kappa_j_part_sum( double p, double e, double J_p, double J_p_prime )
{
  double kappa_j = 0;
  double over_e = 1.0/e;
	double e_pow_2 = e*e;
  double over_ee = 1.0/e_pow_2;
  double e_factor = 1.0-e_pow_2;
  double p_pow_2 = p*p;

  kappa_j = 0.5 * p_pow_2 * sqrt(e_factor) * ( 
      ( (-2.0 * over_ee + 3.0) * over_ee - 1.0) * p * J_p*J_p +
      ( 2.0 * (e + (over_ee - 2.0) * over_e ) * p_pow_2 +
        over_e * (-1.0 + 2.0*over_ee) )  * J_p_prime * J_p +
      2.0 * (1.0 - over_ee) * p * J_p_prime * J_p_prime );
  return kappa_j;
}  

double kappa( kappa_type kappa_t, int p_max, double e )
{
  int p;
  double p_times_e;
  double J_p;
  double J_p_next;
  double J_p_prime;

  double result = 0;

  for ( p = 1; p <= p_max; ++p )
  {
    p_times_e = (double)p*e;

    J_p = gsl_sf_bessel_Jn( p, p_times_e );
    J_p_next  = gsl_sf_bessel_Jn( p+1, p_times_e );
    J_p_prime = (p/p_times_e)*J_p - J_p_next;

    if ( kappa_t == kappa_e )
    {
      result += kappa_e_part_sum( (double) p, e, J_p, J_p_prime );
    }
    else if ( kappa_t == kappa_j )
    {
      result += kappa_j_part_sum( (double) p, e, J_p, J_p_prime );
    }
    else
    {
      fprintf( stderr, "Error in kappa_t: %d\n", kappa_t);
      exit(1);
    }
  }
  return result;
}

double x_dot_0pn( double e, double eta ) /* Eq. (A26) */
{
  double e_pow_2 = e*e;
  double e_factor = 1.0 - e_pow_2;
  return ( 2.* eta * ( (37.*e_pow_2 + 292.)*e_pow_2 + 96. ) ) / ( 15.* pow(e_factor, 7./2.) );
}

double x_dot_1pn( double e, double eta ) /* Eq. (A27) */
{
  double e_factor = 1. - e*e;
  return ( (eta / (420.* pow(e_factor, 9./2.))) *
      ( -(8288.* eta - 11717.) * pow(e, 6.)
        - 14.* (10122.* eta - 12217.) * pow(e, 4.)
        - 120.* (1330.* eta - 731.) * e*e 
        - 16.* (924.* eta + 743.) ) );
}

double x_dot_1p5pn( double e, double eta ) /* Eq. (A28) */
{
  const int p_max = 5; /* upper limit of summation */
  double xdot_1p5;

  if ( e )
  {
    xdot_1p5 = (256.0/5.0) * eta * M_PI * kappa( kappa_e, p_max, e );
  }
  else 
  {
    xdot_1p5 = (256.0 * M_PI * eta / 5.0 ); /* T4 1.5pN term */
  }
   
  return xdot_1p5;
}

double x_dot_2pn( double e, double eta ) /* Eq. (A29) */
{
  double e_factor = 1. - e*e;
  double e_pow_2 = e*e;
  double e_pow_4 = e*e*e*e;
  double eta_pow_2 = eta*eta;

  return ( eta / (45360.* pow(e_factor, 11./2.)) ) * ( 
      ( (1964256.* eta_pow_2 - 3259980.* eta + 3523113.) * e_pow_4 + 
        (64828848.* eta_pow_2 - 123108426.* eta + 83424402.)  * e_pow_2 + 
        (16650606060.* eta_pow_2 - 20720426.* eta + 783768.) ) * e_pow_4 +
        (61282032.* eta_pow_2 + 15464736.* eta - 92846560.) * e_pow_2 +
          1909104.* eta + sqrt(e_factor) * ( ( (264600. - 1058400.* eta) * e_pow_4 +
          (64532160. - 25812864.* eta) ) * e_pow_2 - 580608.*eta + 1451520.) +
        4514976.* eta - 360224. );
}
#if 0

/* Full 3.5pN T4 */
double x_dot_t4( double x )
{
  double xdot;
  double x_pow_5 = x*x*x*x*x;
  double x_factor = (16.0/5.0) * x_pow_5;
  const double euler_gamma = 0.57721566490153286;

 xdot = (16.0 * x_pow_5 / 5.0) * (1.0 - ( 487./168.0 +
        4.0*M_PI * sqrt(x) ) * x + ( 274229.0/72576.0 -    
        254.*M_PI*sqrt(x) ) * x*x + ( 178384023737.0 / 3353011200.0 + 
        1475.0 * M_PI*M_PI / 192.0 - 1712.0*euler_gamma / 105.0 -
        856.0 * ln(16.0*x) + 3310.0 * M_PI * sqrt(x) / 189.0 ) * x*x*x;
  return xdot;
}
#endif


double e_dot_0pn( double e, double eta ) /* Eq. (A31) */
{
	double e_pow_2 = e*e;
  double e_factor = 1. - e_pow_2;
  return  -( e * eta * ( 121.* e_pow_2 + 304. )) /
      ( 15.* pow(e_factor, 5./2.));
}

double e_dot_1pn( double e, double eta ) /* Eq. (A32) */
{
  double e_pow_2 = e*e;
  double e_factor = 1.0 - e_pow_2;

  return ((e*eta) / (2520.* pow(e_factor, 7./2.))) * 
      ( ( (93184.* eta - 125361.) * e_pow_2 +
          12.* (54271.*eta - 59834.) ) * e_pow_2 +
        8.* (2588.*eta + 8451.) );
}

double e_dot_1p5pn( double e, double eta ) /* Eq. (A33) */
{
  const int p_max = 5;
  double e_pow_2 = e*e;
  double e_factor = 1.0 - e_pow_2;
  double edot_1p5;

  if ( e )
  {
    edot_1p5 =  ( (128.0 * eta * M_PI)/(5.0*e) ) * 
      ( (e_pow_2 - 1.0) * kappa(kappa_e, p_max, e) +
        sqrt(e_factor) * kappa(kappa_j, p_max, e) );
  }  
  else 
  {  
    edot_1p5 = 0.0;
  }

  return edot_1p5;
}

double e_dot_2pn( double e, double eta ) /* Eq. (A34) */
{
  double eta_pow_2 = eta*eta;
  double e_pow_2 = e*e;
  double e_pow_4 = e_pow_2*e_pow_2;
  double e_factor = 1. - e_pow_2;

  return -( (e*eta) / (30240.* pow(e_factor, 9./2.)) ) *
      ( ( (2758560.* eta_pow_2 - 4344852.* eta + 3786543.) * e_pow_4 +
          (42810096.* eta_pow_2 - 78112266.* eta + 46579718.) * e_pow_2 +
        (48711348.* eta_pow_2  - 35583228.* eta - 36993396.) ) * e_pow_2 +
        4548096.* eta_pow_2 + sqrt(e_factor) *
        ( ( (2847600. - 1139040.* eta) * e_pow_2 + 
            (35093520. - 14037408.* eta) ) * e_pow_2 - 
          5386752.* eta  + 13466880. ) +
            13509360. *eta - 15198032. );
}

double l_dot_1pn( double e, double eta ) /* Eq. (A2) */
{
  return ( 3. / (e*e - 1.) );
}

double l_dot_2pn( double e, double eta ) /* Eq. (A3) */
{
  return ( ((26.*eta - 51.)*e*e + 28. - 18.) /
      (4.*(e*e-1.) * (e*e-1.)) );
}

double l_dot_3pn( double e, double eta ) /* Eq. (A4) */
{

  return ( (-1./(128.*pow(1.-e*e,7./2.))) * 
      ((1536.*eta*eta*eta-3840.) *e*e*e*e + (1920.-786.*eta) *e*e
       -768.*eta + sqrt(1.-e*e) * 
       ( (1040.*eta*eta-1760.*eta+2496.)*e*e*e*e
         + (5120.*eta*eta+123.*eta*M_PI*M_PI-17856.*eta+8544)*e*e
         +896.*eta-14624.*eta+492.*eta*M_PI*M_PI-192.) + 1920.));
}

static double cosu_factor( double e, double u ) 
{
  return (e * cos(u) - 1.);
}

double phi_dot_0pn( double e, double eta, double u ) /* Eq. (A11) */
{
  return sqrt( 1.0 - e*e ) / ( cosu_factor(e, u) * cosu_factor(e, u) );
}                  	

double phi_dot_1pn( double e, double eta, double u ) /* Eq. (A12) */
{
  double u_factor = cosu_factor(e, u);

  return -( e * ( eta - 4.0 ) * ( e - cos(u) ) /
      ( sqrt( 1.0 - e*e ) * u_factor*u_factor*u_factor ) );
}

double phi_dot_2pn( double e, double eta, double u ) /* Eq. (A13) */
{
  double u_factor = cosu_factor(e, u);
  double e_pow_2 = e*e;
  double e_pow_4 = e_pow_2 * e*e;
  double e_pow_6 = e_pow_4 * e*e;
  double e_factor = 1.0 - e_pow_2;
  double eta_pow_2 = eta*eta;
  double cos_u_pow_2 = cos(u) * cos(u);

  return ( (1. / (12.* pow(e_factor, 3./2.) * pow(u_factor, 5.))) *
      ( (-12.* eta - 18.) * eta * e_pow_6 + 
        (20. * eta_pow_2 - 26. * eta - 60.) * e_pow_4 +
        ( -2.* eta_pow_2 + 50.* eta + 75. ) * e_pow_2 + 
        ( ( -14.* eta_pow_2 + 8.* eta - 147.) * e*e_pow_4 +
          (8. * eta_pow_2 + 22.* eta + 42. )* e * e_pow_2 ) * cos(u)*cos_u_pow_2 + 
        ((17.* eta_pow_2 - 17.* eta + 48.) * e_pow_6 + 
         ( -4.* eta_pow_2 -38.* eta + 153. ) * e_pow_4 +
         (5.* eta_pow_2 - 35.* eta + 114.)* e_pow_2 ) * cos_u_pow_2 - 36.* eta + 
        ( ( -eta_pow_2 + 97.* eta + 12.) * e_pow_4*e + 
         ( -16.* eta_pow_2 - 74.* eta - 81.) * e * e_pow_2 +
         ( -eta_pow_2 + 67.* eta - 246.) * e ) * cos(u) + 
        sqrt( e_factor ) * ( e * e_pow_2 * (36.* eta - 90.) * cos(u)*cos_u_pow_2 +
          ( ( 180. - 72.* eta) * e_pow_4 + 
            ( 90. - 36.* eta) * e_pow_2 ) * cos_u_pow_2 + 
          ( ( 144.* eta - 360.) * e*e_pow_2 + 
            ( 90. - 36.* eta ) * e ) * cos(u) +
          ( 180. - 72. * eta ) * e_pow_2 + 36.* eta - 90.) + 90.) );
}

double phi_dot_3pn( double e, double eta, double u )
{
  double u_factor = cosu_factor(e, u);
  double pi_pow_2 = M_PI*M_PI;
  double eta_pow_2 = eta*eta;
  double e_pow_2 = e*e;
  double e_pow_4 = e*e * e_pow_2;
  double e_factor = 1.0-e_pow_2;

  return ( (1./ ( 13440.*pow(e_factor, 3./2.) * pow( u_factor, 7.) )) *
      ( ( (10080.*eta_pow_2 + 40320.*eta - 15120.) * eta * e_pow_4*e_pow_4 +
        (-52640.*eta_pow_2 - 13440.*eta + 483280.) * eta * e_pow_4*e_pow_2 +
        ( (84000.*eta_pow_2 - 190400.*eta - 17220.*pi_pow_2 - 50048.) * eta - 241920.) * e_pow_4 +
        ( (-52640.*eta_pow_2 + 516880.*eta + 68880.*pi_pow_2 - 1916048.) * eta + 262080.) * e_pow_2 +
        ( (4480.*eta_pow_2 - 412160.*eta  - 30135.*pi_pow_2 + 553008.) * eta + 342720.) ) * e_pow_2 + 
        ( ( ( (13440.*eta_pow_2 + 94640.*eta - 113680.) * eta - 221760.) * e_pow_4 +
         ( (-11200.*eta_pow_2 - 112000.*eta + 12915.*pi_pow_2 + 692928.) * eta - 194880.) * e_pow_2 + 
         ( (4480.*eta_pow_2 + 8960.*eta - 43050.*pi_pow_2 + 1127280.) * eta - 147840.) ) * e * e_pow_4 ) * pow( cos(u), 5.) +
        ( ( (-16240.*eta_pow_2 + 12880.*eta + 18480.) * eta * e_pow_4 * e_pow_2 +
         ( (16240.*eta_pow_2 - 91840.*eta + 17220.*pi_pow_2 - 652192.) * eta + 100800.) * e_pow_4 +
         ( (-55440.*eta_pow_2 + 34160.*eta - 30135.*pi_pow_2 - 2185040.) * eta + 2493120.) * e_pow_2 + 
         ( (21480.*eta_pow_2 + 86800.*eta + 163590.*pi_pow_2 - 5713888.) * eta + 228480.) ) * e_pow_4 ) * pow(cos(u), 4.) + 
        ( ( ( (560.*eta_pow_2 - 137200.*eta + 388640.) * eta + 241920.) * e_pow_4 * e_pow_2 + 
          ( (30800.*eta_pow_2 - 264880.*eta - 68880.*pi_pow_2 + 624128.) * eta + 766080.) * e_pow_4 + 
          ( (66640.*eta_pow_2 + 612080.*eta - 8610.*pi_pow_2 + 6666080. ) * eta - 6652800.) * e_pow_2 + 
         ( (-30800.*eta_pow_2-364880.*eta- 223860.*pi_pow_2 + 9386432.) * eta) ) * e*e_pow_2 ) * pow(cos(u), 3.) +
        67200. * eta_pow_2 + 
        ( ( (4480.*eta_pow_2 - 20160.*eta + 16800.) * eta * e_pow_4*e_pow_4 +
       ( (3920.*eta_pow_2 + 475440.*eta - 17220.*pi_pow_2 + 831952.) * eta - 7257600.) * e_pow_4*e_pow_2 + 
            ( (-75600.*eta_pow_2 + 96880.*eta + 154980.*pi_pow_2 - 3249488.) * eta - 685440.) * e_pow_4 +
            ( (5040.*eta_pow_2 - 659120.*eta + 25830.*pi_pow_2 - 7356624.) * eta + 6948480.) * e_pow_2 +
            ( (-5040.*eta_pow_2 + 190960.*eta + 137760.*pi_pow_2 - 7307920.) * eta + 107520.) ) * e_pow_2 ) * pow(cos(u), 2.) -
        761600.*eta + 
        ( ( ( (-2240.*eta_pow_2 - 168000.*eta - 424480.) * eta * e_pow_4  + 
            ( (28560.*eta_pow_2 + 242480.*eta + 34440.*pi_pow_2 - 1340224.) * eta + 725760.) * e_pow_2 + 
            ( (-33040.* eta_pow_2 - 754880.* eta - 172200.*pi_pow_2 + 5458480. ) * eta - 221760.) ) * e_pow_4 + 
            ( (40880.*eta_pow_2 + 738640.*eta + 30135.*pi_pow_2 + 1554048.) * eta - 2936640.) * e_pow_2 + 
            ( (-560.*eta_pow_2  - 100240.*eta - 43050.*pi_pow_2 + 3284816.) * eta - 389760.) ) * e ) * cos(u) +
        sqrt(e_factor) * ( ( (-127680.*eta*eta + 544320.*eta-739200.) * e_pow_2 +
              ( (-53760.*eta - 8610.*pi_pow_2  + 674240.) * eta - 67200.) * e*e_pow_4 ) * pow(cos(u), 5.) + 
            ( ( ( (161280.*eta - 477120.) * eta + 537600.) * e_pow_4 + 
              ( (477120.*eta + 17220.*pi_pow_2 - 2894080.) * eta + 2217600.) * e_pow_2 +
               ( (268800.*eta + 25830.* pi_pow_2 - 2721600.) * eta + 1276800.) ) * e_pow_4 ) * pow(cos(u),4.) + 
            (( (-524160.*eta + 1122240.)*eta - 940800.) * pow(e,7.) + 
             ( (-873600.*eta - 68880.* pi_pow_2 + 7705600.) * eta - 3897600.) * e*e*e*e*e + 
             ( (-416640.*eta - 17220.* pi_pow_2 + 3357760.) * eta - 
                 3225600.) * e*e*e ) * cos(u)*cos(u)*cos(u) + 
            ((604800.*eta*eta - 50400.*eta - 403200.) * pow(e,6.) 
             + ( (1034880.*eta + 103320.*pi_pow_2 - 11195520.) * eta + 5779200.) * e*e*e*e +
             ( (174720.*eta - 17229.*pi_pow_2 - 486080.) * eta + 2688000.) * e*e ) * cos(u)*cos(u)
            + ((-282240.*eta_pow_2 - 450240.*eta + 1478400.) *e*e*e*e*e
              + ( (-719040.*eta - 68880*pi_pow_2 + 8128960.) * eta - 5040000.) * e*e*e +
              ( (94080.*eta + 25830.*pi_pow_2 - 1585920.) * eta - 470400.) * e ) * cos(u) 
            - (67200.*eta + 761600.)*eta + e*e*e*e * (40320.*eta*eta + 309120.*eta -672000.) + e*e * (208320.*eta*eta + 
                17220.* pi_pow_2*eta - 2289280.*eta + 
                1680000.) - 8610.*pi_pow_2*eta - 201600.) + 8610.*pi_pow_2*eta + 201600.) );
}
