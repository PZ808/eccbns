/* $Id: x_ecc_parse_args.c,v 1.1 2009/06/02 18:14:19 pzimmerman Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <dlfcn.h>
#include "x_ecc_parse_args.h"
#include "x_drive.h"

extern int con_pn_order;
extern int rad_pn_order;
extern double mass1;
extern double mass2;
extern double f_gw_init;
extern double e_init;
extern double mean_anom_init;
extern double ode_eps;
extern double sampling_rate;

int parse_args( int argc, char *argv[] )
{
  struct option long_options[] =
  {
    {"cPN-order", required_argument,  0,  'C'},
    {"rPN-order", required_argument,  0,  'R'},
    {"mass-1",	  required_argument,	0,	'm'},
    {"mass-2",	  required_argument,	0,	'M'},
    {"f-init",	  required_argument,	0,	'f'},
    {"e-init",  	required_argument,	0,	'e'},
    {"l-init",  	required_argument,	0,	'l'},
    {"accuracy",	required_argument,	0,	'a'},
    {"samp-rate",	required_argument,	0,	'r'},
    {0, 0, 0, 0}
  };

  int c;

  while( 1 )
  {
    int option_index = 0;
    c = getopt_long_only( argc, argv, 
        "C:R:m:M:f:e:l:a:r", long_options, &option_index );

    if ( c == -1 ) break;

    switch ( c )
    {
      case 0:
        /* if this option set a flag, do nothing else now */
        if ( long_options[option_index].flag != 0 )
        {
          break;
        }
        else
        {
          fprintf( stderr, "error parsing option %s with argument %s\n",
              long_options[option_index].name, optarg );
          exit( 1 );
        }
        break;

      case 'C':
        con_pn_order = (int) atoi( optarg );
        break;

      case 'R':
        rad_pn_order = (int) atoi( optarg );
        break;

      case 'm':
        mass1 = (double) atof( optarg );
        break;

      case 'M':
        mass2 = (double) atof( optarg );
        break;

      case 'f':
        f_gw_init = (double) atof( optarg );
        break;

      case 'e':
        e_init = (double) atof( optarg );
        break;

      case 'l':
        mean_anom_init = (double) atof( optarg );
        break;

      case 'a':
        ode_eps = (double) atof( optarg );
        break;

      case 'r':
        sampling_rate = (double) atof( optarg );
        break;

      default:
        fprintf( stderr, "Error: unknown error while parsing options (%d)\n", c );
        exit( 1 );
    }
  }


  if ( optind < argc )
  {
    fprintf( stderr, "Error: extraneous command line arguments:\n" );
    while ( optind < argc )
    {
      fprintf ( stderr, "%s\n", argv[optind++] );
    }
    exit( 1 );
  }

  if ( (con_pn_order < 0) || (con_pn_order > 3) )
  {
    fprintf( stderr, "Error: invalid argument to --cPN-order (%d)\n", con_pn_order );
    exit( 1 );
  }

  if ( (rad_pn_order < 0) || (rad_pn_order > 2) )
  {
    fprintf( stderr, "Error: invalid argument to --rPN-order (%d)\n", rad_pn_order );
    exit( 1 );
  }

  if ( mass1 < 0 )
  {
    fprintf( stderr, "Error: invalid argument to --mass-1 (%f)\n", mass1 );
    exit( 1 );
  }

  if ( mass2 < 0 )
  {
    fprintf( stderr, "Error: invalid argument to --mass-2 (%f)\n", mass2 );
    exit( 1 );
  }

  if ( (ode_eps <= 0) || (ode_eps < 1.0e-16) )
  {
    fprintf( stderr, "Error: invalid argument to --accuracy (%f)\n", ode_eps );
    exit( 1 );
  }

  if ( e_init < 0 || e_init > 1)
  {
    fprintf( stderr, "Error: invalid argument to --e-init (%f)\n", e_init );
    exit( 1 );
  }

  if ( f_gw_init <= 0 )
  {
    fprintf( stderr, "Error: invalid argument to --f-init (%f)\n", f_gw_init );
    exit( 1 );
  }

  if ( mean_anom_init < 0 )
  {
    fprintf( stderr, "Error: invalid argument to --l-init (%f)\n", mean_anom_init );
    exit( 1 );
  }

  return 0;
}
