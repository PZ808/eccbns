# Makefile for x-based pn codes
#
# $Id: Makefile,v 1.12 2009/07/16 03:04:35 pzimmerman Exp $
#

CC = gcc
CFLAG = -Wall 
CPATH = -I. -I/opt/lscsoft/non-lsc/include
LIBS = -lgsl -lgslcblas 
LPATH = -L/opt/lscsoft/non-lsc/lib

OBJA = x_ecc_kepler.o x_ecc_ode_pn_terms.o x_ecc_ode.o x_ecc_sep.o 
SRCA = x_ecc_kepler.c x_ecc_ode_pn_terms.c x_ecc_ode.c x_ecc_sep.c 

OBJB = x_ecc_kep_newt.o x_ecc_ode_pn_terms.o x_ecc_ode_newt.o x_ecc_sep.o 
SRCB = x_ecc_kep_newt.c x_ecc_ode_pn_terms.c x_ecc_ode_newt.c x_ecc_sep.c 

all: x_drive NumCycles xwave_dev

x_ecc_wave: x_ecc_wave.o $(OBJA) 
	$(CC) $(LPATH) -pg -o x_ecc_wave x_ecc_wave.o $(OBJA) $(LIBS) -lm 

x_ecc_wave.o: x_ecc_wave.c $(SRCA) 
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_wave.c $(SRCA)

x_drive: x_drive.o $(OBJB) 
	$(CC) $(LPATH) -pg -o x_drive x_drive.o $(OBJB) $(LIBS) -lm

x_drive.o: x_drive.c $(SRCB)
	$(CC) $(CFLAG)$(CPATH) -pg -c x_drive.c $(SCRB)

xwave_dev: xwave_dev.o $(OBJB) 
	$(CC) $(LPATH) -pg -o ecc_hoft  xwave_dev.o $(OBJB) $(LIBS) -lm

xwave_dev.o: xwave_dev.c $(SRCB)
	$(CC) $(CFLAG)$(CPATH) -pg -c xwave_dev.c $(SCRB)

NumCycles: NumCycles.o $(OBJB)
	$(CC) $(LPATH) -pg -o NumCycles NumCycles.o $(OBJB) $(LIBS) -lm

NumCycles.o: NumCycles.c
	$(CC) $(CFLAG)$(CPATH) -pg -c NumCycles.c $(SRCB)

x_ecc_ode.o: x_ecc_ode.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_ode.c

x_ecc_ode_newt.o: x_ecc_ode_newt.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_ode_newt.c

x_ecc_ode_pn_terms.o: x_ecc_ode_pn_terms.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_ode_pn_terms.c

x_ecc_sep.o: x_ecc_sep.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_sep.c

x_ecc_kepler.o: x_ecc_kepler.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_kepler.c

x_ecc_kep_newt.o: x_ecc_kep_newt.c
	$(CC) $(CFLAG)$(CPATH) -pg -c x_ecc_kep_newt.c

clean:
	rm -f x_ecc_wave x_ecc_wave.o x_drive x_drive.o NumCycles NumCycles.o $(OBJA) $(OBJB)

clean-data:
	rm -f xwave_[0-3]_[0-4]_*.txt xdyn_[0-3]_[0-4]_*.txt xtraj_[0-3]_[0-4]_*.txt


