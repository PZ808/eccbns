#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <dlfcn.h>
#include "parse_args.h"
#include "ode_system.h"

extern double t_min;
extern double t_max;
extern double t_step;
extern double t_scale;
extern double eps;
extern double hguess;
extern char solver_type[];
extern char ode_system_name[];
extern double e_init;
extern double f_init;
extern double m1;
extern double m2;

extern int (*init_y)(double, double[], void* );
extern int (*ode_system)(double, const double[], double dydt[], void* );
extern int (*ode_term_cond)(double, double y[], void* );
extern int (*check_args)(void);
extern int (*ode_filename)(void);
extern int (*write_output)( FILE*, double, double, double[], void* );

int parse_args( int argc, char *argv[] )
{
	void *ode_system_library;
	char ode_system_library_name[ODE_SYSTEM_NAME_MAX];
	int have_solver_type = 0;
  int have_ode_system_name = 0;
  int have_hguess = 0;
	struct option long_options[] =
	{
		{"t-min",	required_argument,	0,	't'},
		{"t-max",	required_argument,	0,	'T'},
		{"t-step",	required_argument,	0,	's'},
		{"t-scale",	required_argument,	0,	'c'},
		{"accuracy",	required_argument,	0,	'a'},
		{"h-guess",	required_argument,	0,	'h'},
		{"solver-type",	required_argument,	0,	'S'},
		{"ode-system",	required_argument,	0,	'o'},
		{"e-init",  	required_argument,	0,	'e'},
		{"f-init",	required_argument,	0,	'f'},
		{"mass-1",	required_argument,	0,	'm'},
		{"mass-2",	required_argument,	0,	'M'},
		{0, 0, 0, 0}
	};

	int c;

	while( 1 )
	{
		int option_index = 0;
		c = getopt_long_only( argc, argv, 
				"t:T:s:c:a:h:S:o:e:f:m:M", long_options, &option_index );

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

			case 't':
				t_min = (double) atof( optarg );
				break;

			case 'T':
				t_max = (double) atof( optarg );
				break;

			case 's':
				t_step = (double) atof( optarg );
				break;

			case 'c':
				t_scale = (double) atof( optarg );
				break;

			case 'a':
				eps = (double) atof( optarg );
				break;

			case 'h':
				hguess = (double) atof( optarg );
                                have_hguess = 1;
				break;

			case 'S':
				strncpy( solver_type, optarg, SOLVER_TYPE_MAX * sizeof(char) );
                                have_solver_type = 1;
				break;

			case 'o':
				strncpy( ode_system_name, optarg, ODE_SYSTEM_NAME_MAX * sizeof(char) );
                                snprintf( ode_system_library_name, ODE_SYSTEM_NAME_MAX * sizeof(char), "./lib%s.so", optarg );
                                have_ode_system_name = 1;
				break;

			case 'e':
				e_init = (double) atof( optarg );
				break;

			case 'f':
				f_init = (double) atof( optarg );
				break;

			case 'm':
				m1 = (double) atof( optarg );
				break;

			case 'M':
				m2 = (double) atof( optarg );
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

	if ( t_min < 0 )
	{
		fprintf( stderr, "Error: invalid argument to --t-min (%f)\n", t_min );
		exit( 1 );
	}

	if ( t_max < 0 )
	{
		fprintf( stderr, "Error: invalid argument to --t-max (%f)\n", t_max );
		exit( 1 );
	}

	if ( t_scale <= 0 )
	{
		fprintf( stderr, "Error: invalid argument to --t-scale (%f)\n", t_scale );
		exit( 1 );
	}

	if ( t_max < t_min )
	{
		fprintf( stderr, "Error: t-max is less than t-min\n" );
		exit( 1 );
	}

	if ( t_step <= 0 )
	{
		fprintf( stderr, "Error: invalid argument to --t-step (%f)\n", t_step );
		exit( 1 );
	}

	if ( eps <= 0 )
	{
		fprintf( stderr, "Error: invalid argument to --accuracy (%f)\n", eps );
		exit( 1 );
	}

        if ( have_hguess )
        {
          if ( hguess <= 0 )
          {
            fprintf( stderr, "Error: invalid argument to --h-guess (%f)\n", hguess );
            exit( 1 );
          }
        }

	if ( ! have_solver_type )
	{
		fprintf( stderr, "Error: --solver-type must be specified\n" );
		exit( 1 );
	}

	if ( ! have_ode_system_name )
	{
		fprintf( stderr, "Error: --ode-system must be specified\n" );
		exit( 1 );
	}
	
        /* try and load the ode system */
        ode_system_library = dlopen( ode_system_library_name, RTLD_LAZY );
        if ( ! ode_system_library )
        {
		fprintf( stderr, "Error opening library: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

        init_y = dlsym( ode_system_library, "init_y" );
        if ( !init_y )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	ode_system = dlsym( ode_system_library, "ode_system" );
        if ( !ode_system )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	write_output = dlsym( ode_system_library, "write_output" );
        if ( !write_output )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	ode_term_cond = dlsym( ode_system_library, "ode_term_cond" );
        if ( !ode_term_cond )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	check_args = dlsym( ode_system_library, "check_args" );
        if ( !check_args )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	ode_filename = dlsym( ode_system_library, "ode_filename" );
        if ( !ode_filename )
        {
		fprintf( stderr, "Error opening function: %s\n", dlerror() );
        	exit( 1 );
        }
        dlerror();

	if ( check_args() ) exit( 1 );

	return 0;
}
