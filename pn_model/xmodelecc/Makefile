# Makefile for x-based pn codes
#
# $Id: Makefile,v 1.2 2009/07/16 02:59:43 pzimmerman Exp $
#

CC = gcc
CFLAG = -Wall 
CPATH = -I. -I/opt/lscsoft/non-lsc/include
LIBS = -lgsl -lgslcblas 
LPATH = -L/opt/lscsoft/non-lsc/lib

OBJ = xecc_kep.o xecc_pn_terms.o xecc_ode.o xecc_sep.o 
SRC = xecc_kep.c xecc_pn_terms.c xecc_ode.c xecc_sep.c 

all: xecc

xecc: xecc.o $(OBJ) 
	$(CC) $(LPATH) -pg -o xecc xecc.o $(OBJ) $(LIBS) -lm

xecc.o: xecc.c $(SRC)
	$(CC) $(CFLAG)$(CPATH) -pg -c xecc.c $(SCRB)

xecc_ode.o: xecc_ode.c
	$(CC) $(CFLAG)$(CPATH) -pg -c xecc_ode.c

xecc_pn_terms.o: xecc_pn_terms.c
	$(CC) $(CFLAG)$(CPATH) -pg -c xecc_pn_terms.c

xecc_sep.o: xecc_sep.c
	$(CC) $(CFLAG)$(CPATH) -pg -c xecc_sep.c

xecc_kepler.o: xecc_kep.c
	$(CC) $(CFLAG)$(CPATH) -pg -c xecc_kep.c

clean:
	rm -f xecc xecc.o $(OBJ) 

clean-data:


