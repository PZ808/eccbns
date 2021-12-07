#####################################################
#                                                   #
# Makefile for orbit.c, wave.c                      #
#                                                   #
#####################################################

# Compiler, options, and libraries
CC = gcc
OPTIONS= -g -Wall -I/home/dbrown/local/nr/include
LIBS= -lrecipes_c
LOPTIONS= -L/home/dbrown/local/nr/lib64

# Optimization
OPTIM= -O6

# Make
all: duration orbit wave fourier ambiguity ff

duration: duration.o
	$(CC) $(OPTIM) -o duration duration.o $(LOPTIONS) $(LIBS) -lm  
duration.o: duration.c
	$(CC) $(OPTIM) $(OPTIONS) -c duration.c

orbit: orbit.o
	$(CC) $(OPTIM) -o orbit orbit.o $(LOPTIONS) $(LIBS) -lm 
orbit.o: orbit.c
	$(CC) $(OPTIM) $(OPTIONS) -c orbit.c

wave: wave.o
	$(CC) $(OPTIM) -o wave wave.o $(LOPTIONS) $(LIBS) -lm  
wave.o: wave.c
	$(CC) $(OPTIM) $(OPTIONS) -c wave.c

fourier: fourier.o
	$(CC) $(OPTIM) -o fourier fourier.o $(LOPTIONS) $(LIBS) -lm  
fourier.o: fourier.c
	$(CC) $(OPTIM) $(OPTIONS) -c fourier.c

ambiguity: ambiguity.o
	$(CC) $(OPTIM) -o ambiguity ambiguity.o $(LOPTIONS) $(LIBS) -lm  
ambiguity.o: ambiguity.c
	$(CC) $(OPTIM) $(OPTIONS) -c ambiguity.c

ff: ff.o
	$(CC) $(OPTIM) -o ff ff.o $(LOPTIONS) $(LIBS) -lm  
ff.o: ff.c
	$(CC) $(OPTIM) $(OPTIONS) -c ff.c

clean: 
	rm -f duration duration.o orbit orbit.o 
	rm -f wave wave.o fourier fourier.o   
	rm -f ambiguity ambiguity.o ff ff.o


