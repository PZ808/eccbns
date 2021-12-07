all: pn_evolve x_formalism_sim x_formalism_new
           
pn_merc: pn_merc.o conservative_dynamics.o
	condor_compile gcc -L/home/pzimmerman/lib -o pn_merc pn_merc.o mikkola.o conservative_dynamics.o -lgsl -lgslcblas -lm

pn_evolve: pn_evolve.o conservative_dynamics.o
	condor_compile gcc -L/home/pzimmerman/lib -o pn_evolve pn_evolve.o mikkola.o conservative_dynamics.o -lgsl -lgslcblas -lm

x_con: x_con.o conservative_dynamics.o conservative_ode_system.o
	gcc -L/home/pzimmerman/lib -o run_x_con x_con.o mikkola.o conservative_dynamics.o conservative_ode_system.o -lgsl -lgslcblas -lm

x_formalism: x_formalism.o conservative_dynamics.o reactive_dynamics.o mean_motion_init.o
	condor_compile gcc -L/home/pzimmerman/lib -o run_x_model x_formalism.o mean_motion_init.o mikkola.o conservative_dynamics.o reactive_dynamics.o -lgsl -lgslcblas -lm

x_formalism_sim: x_formalism_sim.o conservative_dynamics.o reactive_dynamics.o 
	condor_compile gcc -L/home/pzimmerman/lib -o run_x_model_sim x_formalism_sim.o mikkola.o conservative_dynamics.o reactive_dynamics.o -lgsl -lgslcblas -lm

x_formalism_new: x_formalism_new.o conservative_dynamics.o reactive_dynamics.o mean_motion_init.o
	condor_compile gcc -L/home/pzimmerman/lib -o run_x_model_new x_formalism_new.o mean_motion_init.o mikkola.o conservative_dynamics.o reactive_dynamics.o -lgsl -lgslcblas -lm

conservative_dynamics.o: conservative_dynamics.c 
	gcc -Wall -I/home/pzimmerman/include -c conservative_dynamics.c

reactive_dynamics.o: reactive_dynamics.c 
	gcc -Wall -I/home/pzimmerman/include -c reactive_dynamics.c

pn_merc.o: pn_merc.c mikkola.c
	gcc -Wall -I/home/pzimmerman/include -c pn_merc.c mikkola.c conservative_dynamics.c

pn_evolve.o: pn_evolve.c mikkola.c
	gcc -Wall -I/home/pzimmerman/include -c pn_evolve.c mikkola.c conservative_dynamics.c

x_con.o: x_con.c conservative_dynamics.c mikkola.o
	gcc -Wall -I/home/pzimmerman/include -c x_con.c mikkola.c conservative_dynamics.c conservative_ode_system.c

x_formalism.o: x_formalism.c mikkola.c conservative_dynamics.c reactive_dynamics.c mean_motion_init.c
	gcc -Wall -I/home/pzimmerman/include -c x_formalism.c mikkola.c conservative_dynamics.c reactive_dynamics.c mean_motion_init.c

x_formalism_sim.o: x_formalism_sim.c mikkola.c conservative_dynamics.c reactive_dynamics.c 
	gcc -Wall -I/home/pzimmerman/include -c x_formalism_sim.c mikkola.c conservative_dynamics.c reactive_dynamics.c

x_formalism_new.o: x_formalism_new.c mikkola.c conservative_dynamics.c reactive_dynamics.c mean_motion_init.c
	gcc -Wall -I/home/pzimmerman/include -c x_formalism_new.c mikkola.c conservative_dynamics.c reactive_dynamics.c mean_motion_init.c

mikkola.o: mikkola.c
	gcc -L/home/pzimmerman/include -c mikkola.c 

mean_motion_init.o: mean_motion_init.c
	gcc -L/home/pzimmerman/include -c mean_motion_init.c 

clean: 
	rm pn_merc pn_merc.o conservatice_dynamics.o mikkola.o pn_evolve pn_evolve.o x_formalism.o run_x_model x_con.o run_x_con run_x_model_new run_x_formalism_simp.o run_x_formalism_simp
