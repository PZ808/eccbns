# Makefile for test ODE systems
#
# $Id: Makefile,v 1.12 2008/07/01 21:04:34 dbrown Exp $
#

all: gsl_solver nr_solver nr_solver_double analytic_soln libtoy.so libeccbns.so libeccbns_physical.so

gsl_solver: gsl_solver.o parse_args.o
	gcc -rdynamic -g -o gsl_solver gsl_solver.o parse_args.o -lgsl -lgslcblas -lm -ldl

nr_solver: nr_solver.o parse_args.o
	gcc -rdynamic -g -O6 -o nr_solver nr_solver.o parse_args.o -L/home/dbrown/local/nr/lib64 -lrecipes_c -lm -ldl

nr_solver_double: nr_solver_double.o parse_args.o nr_double.o
	gcc -rdynamic -g -O6 -o nr_solver_double nr_solver_double.o nr_double.o parse_args.o -lm -ldl

analytic_soln: analytic_soln.o
	gcc -g -o analytic_soln analytic_soln.o -lm

libeccbns.so: eccbns_ode_system.o
	gcc -g -fPIC -nostartfiles -shared -o libeccbns.so eccbns_ode_system.o -lm

eccbns_ode_system.o: eccbns_ode_system.c
	gcc -g -fPIC -Wall -o eccbns_ode_system.o -c eccbns_ode_system.c

libeccbns_physical.so: eccbns_physical_ode_system.o
	gcc -g -fPIC -nostartfiles -shared -o libeccbns_physical.so eccbns_physical_ode_system.o -lm

eccbns_physical_ode_system.o: eccbns_physical_ode_system.c
	gcc -g -fPIC -Wall -o eccbns_physical_ode_system.o -c eccbns_physical_ode_system.c

libtoy.so: toy_ode_system.o
	gcc -g -fPIC -nostartfiles -shared -o libtoy.so toy_ode_system.o -lm

toy_ode_system.o: toy_ode_system.c
	gcc -g -fPIC -Wall -o toy_ode_system.o -c toy_ode_system.c

parse_args.o: parse_args.c
	gcc -g -Wall -c parse_args.c

gsl_solver.o: gsl_solver.c
	gcc -g -Wall -c gsl_solver.c

nr_solver.o: nr_solver.c
	gcc -g -O6 -g -Wall -I/home/dbrown/local/nr/include -c nr_solver.c

nr_solver_double.o: nr_solver_double.c
	gcc -g -O6 -g -Wall -c nr_solver_double.c

nr_double.o: nr_double.c
	gcc -g -O6 -g -Wall -c nr_double.c

analytic_soln.o: analytic_soln.c
	gcc -g -Wall -c analytic_soln.c

clean:
	rm -f gsl_solver nr_solver nr_solver_double analytic_soln *.o *.so

clean-data:
	rm -f gsl_solver*.txt nr_solver*.txt analytic_soln.txt
