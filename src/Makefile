all: evolve_25pn n_formalism evolve2Rad2ConPN eT_35pn

evolve_25pn: evolve_25pn.o mean_motion_init.o
	gcc -L/opt/lscsoft/non-lsc/lib -o evolve_25pn evolve_25pn.o mean_motion_init.o -lgsl -lgslcblas -lm

evolve2Rad2ConPN: evolve2Rad2ConPN.o mean_motion_init.o
	gcc -L/opt/lscsoft/non-lsc/lib -o n_model_evolution evolve2Rad2ConPN.o mean_motion_init.o -lgsl -lgslcblas -lm
#	condor_complie gcc -L/home/pzimmerman/lib - evolve_n_model_2r2c evolve2Rad2ConPN.o mean_motion_init.o -lgsl -lgslcblas lm 

n_formalism: n_formalism.o mean_motion_init.o
	gcc -L/opt/lscsoft/non-lsc/lib -o run_n_model n_formalism.o mean_motion_init.o -lgsl -lgslcblas -lm

eT_35pn: eT_35pn.o eT_ICs.o
	gcc -L/opt/lscsoft/non-lsc/lib -o eT_35pn eT_35pn.o eT_ICs.o -lgsl -lgslcblas -lm

mean_motion_init.o: mean_motion_init.c
	gcc -Wall -I/opt/lscsoft/non-lsc/include -c mean_motion_init.c

n_formalism.o: n_formalism.c 
	gcc -Wall -I/opt/lscsoft/non-lsc/include -c n_formalism.c 

evolve2Rad2ConPN.o: evolve2Rad2ConPN.c
	gcc -Wall -I/opt/lscsoft/non-lsc/include -c evolve2Rad2ConPN.c 

evolve_25pn.o: evolve_25pn.c
	gcc -Wall -I/opt/lscsoft/non-lsc/include -c evolve_25pn.c

eT_35pn.o: eT_35pn.c eT_ICs.c
	gcc -Wall -I/opt/lscsoft/non-lsc/include -c eT_35pn.c eT_ICs.c
	
clean:
	rm -f evolve2Rad2ConPN evolve2Rad2ConPN.o evolve_25pn evolve_25pn.o  eT_35pn eT_35pn.o
