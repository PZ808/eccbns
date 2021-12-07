/* $Id: x_ecc_sep.c,v 1.4 2009/06/08 18:07:10 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

#include "x_drive.h"

static double cosu_factor(double e, double u)
{
  return (e * cos(u) - 1.);
}

static double rel_sep_0pn( double e, double u )
{
  return ( 1.0 - e*cos(u) );
}

static double rel_sep_1pn( double e, double u, double eta )
{
  double u_factor = cosu_factor(e, u);

  return ( 2.*(u_factor) / (e*e - 1.) +
    (1./6.) * (2.*(eta - 9.) +
    e * (7.*eta - 6.) * cos(u)) );
}

static double rel_sep_2pn(double e, double u, double eta)
{
  double eta_pow_2 = eta*eta;
  double e_pow_2 = e*e;
  double e_factor = 1.0-e*e;

  return (1. / (e_factor*e_factor)) * (1./72.) * (
        ( (8.*eta_pow_2 + 30.*eta + 72.) * e_pow_2 +
          (-16.0*eta_pow_2 - 876.*eta + 756.) ) * e_pow_2 +
        (8.*eta_pow_2 + 198.*eta + 360.) + 
        ( ( ( (-35.*eta_pow_2 + 231.*eta - 72.) * e_pow_2 +
                (70.*eta_pow_2 - 150.*eta - 468.) ) * e_pow_2 +
          (-35.*eta_pow_2 + 567.* eta - 648.) ) * e ) * cos(u) +
        sqrt(e_factor) * ( (360. - 144.*eta)*e_pow_2 + (144.* eta - 360.) +
          ( (180. - 72.*eta) * e_pow_2  + (72.* eta - 180.) )*e*cos(u) ) );
}

static double rel_sep_3pn( double e, double u, double eta )
{
  double pi_pow_2 = M_PI*M_PI;
  double e_factor = 1.0 - e*e;
  double eta_pow_2 = eta*eta;
  double e_pow_2 = e*e;
  double e_pow_4 = e*e*e_pow_2;

  return (1. / (181440.*sqrt(e_factor)*e_factor*e_factor)) *
      ( ( (-665280.*eta_pow_2 + 1753920.* eta - 1814400.) * e_pow_4 +
          ( (725760.*eta_pow_2 - 77490.*pi_pow_2 + 5523840.)*eta - 3628800.) * e_pow_2 +
          ( (544320.*eta_pow_2 + 154980.*pi_pow_2 - 14132160.) * eta + 7257600.) ) * e_pow_2 
            - 604800.*eta_pow_2 + 6854400.*eta +
        ( ( (302400.*eta_pow_2 - 1254960.*eta + 453600.) * e_pow_4 +
            ( (-1542240.*eta_pow_2 - 38745.*pi_pow_2 + 6980400.) * eta - 453600.) * e_pow_2 +
            ( (2177280.*eta_pow_2 + 77490.*pi_pow_2 - 12373200. )*eta + 4989600.) ) * e*e_pow_2 +
          ( (-937440.*eta_pow_2-37845.*pi_pow_2 + 6647760.) * eta - 4989600.) * e ) * cos(u) +
        sqrt(e_factor) * ( ( ( (-4480.*eta_pow_2 - 25200.*eta + 22680.) * eta 
              - 120960.) * e_pow_4 +
            (13440.*eta_pow_2*eta + 4404960.*eta*eta + 116235.* pi_pow_2 - 12718296.* eta + 5261760.) * e_pow_2 +
            ( (-13440.* eta_pow_2 + 2242800.*eta + 348705.* pi_pow_2 - 19225080.) * eta + 1614160.) ) * e_pow_2 +
          (4480.*eta_pow_2 + 45360.*eta - 8600904.) * eta +
          ( ( ( (-6860.*eta_pow_2 - 550620.*eta - 986580.) * eta + 120960.) * e_pow_4 +
              ( (20580.*eta_pow_2 - 2458260.*eta + 3458700.) * eta - 2358720.) * e_pow_2 +
              ( (-20580.*eta_pow_2 - 3539340.*eta - 116235.*pi_pow_2 + 20173860.) * eta - 16148160.) ) * e*e_pow_2 +
          ( (6860.*eta_pow_2 - 1220940.*eta + 464940.*pi_pow_2 + 17875620.) * eta - 417440.) * e ) * cos(u)
          + 116235.*pi_pow_2*eta + 1814400.) - 77490.*pi_pow_2*eta - 1814400. );
}

double separation( int conservative_pn_order,
    double u, double eta, double x, double e) 
{
  if ( conservative_pn_order == 0 )
  {
    return (1.0/x) * rel_sep_0pn( e, u );
  }
  else if  ( conservative_pn_order == 1 )
  {
    return (1.0/x) * rel_sep_0pn( e, u ) + rel_sep_1pn( e, u, eta );
  }
  else if ( conservative_pn_order == 2 )
  {
    return (1.0/x) * rel_sep_0pn( e, u ) + rel_sep_1pn( e, u, eta ) +
      rel_sep_2pn( e, u, eta ) * x;
  }
  else if ( conservative_pn_order == 3 )
  {
    return (1.0/x) * rel_sep_0pn( e, u ) + rel_sep_1pn( e, u, eta ) +
      rel_sep_2pn( e, u, eta ) * x + rel_sep_3pn( e, u, eta ) * x * x;
  }
  else
  {
    fprintf( stderr, "Error in pN order: %d\n", conservative_pn_order );
    exit(1);
  }
}

double x_cartesian( double r, double phi )
{
	return  r * cos(phi);
}

double y_cartesian( double r, double phi )
{
	return  r * sin(phi);
}
