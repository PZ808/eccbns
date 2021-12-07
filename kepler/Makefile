all: eccentric_orbit_newt relative_separation_guide relative_separation_mikkola
           
eccentric_orbit: eccentric_orbit.o kepler.o
	condor_compile gcc -L/home/pzimmerman/lib -o eccentric_orbit eccentric_orbit.o kepler.o -lgsl -lgslcblas -lm
eccentric_orbit_newt: eccentric_orbit_newt.o kepler.o mikkola.o
	condor_compile gcc -L/home/pzimmerman/lib -o eccentric_orbit_newt eccentric_orbit_newt.o kepler.o mikkola.o -lgsl -lgslcblas -lm
pn_ecc_orbit: pn_ecc_orbit.o pn_kepler.o conservative_dynamics.o
	condor_compile gcc -L/home/pzimmerman/lib -o pn_ecc_orbit pn_ecc_orbit.o pn_kepler.o conservative_dynamics.o -lgsl -lgslcblas -lm
relative_separation: relative_separation.o kepler.o mikkola.o
	condor_compile gcc -L/home/pzimmerman/lib -o relative_separation relative_separation.o kepler.o mikkola.o -lgsl -lgslcblas -lm
relative_separation_guide: relative_separation_guide.o kepler.o 
	condor_compile gcc -L/home/pzimmerman/lib -o relative_separation_guide relative_separation_guide.o kepler.o -lgsl -lgslcblas -lm
relative_separation_mikkola: relative_separation_mikkola.o mikkola.o 
	condor_compile gcc -L/home/pzimmerman/lib -o relative_separation_mikkola relative_separation_mikkola.o mikkola.o -lgsl -lgslcblas -lm
kepler.o: kepler.c
	gcc -L/home/pzimmerman/include -c kepler.c 
mikkola.o: mikkola.c
	gcc -L/home/pzimmerman/include -c mikkola.c 
pn_kepler.o: pn_kepler.c
	gcc -L/home/pzimmerman/include -c pn_kepler.c 
conservative_dynamics.o: conservative_dynamics.c 
	gcc -Wall -I/home/pzimmerman/include -c conservative_dynamics.c
eccentric_orbit.o: eccentric_orbit.c
	gcc -Wall -I/home/pzimmerman/include -c eccentric_orbit.c kepler.c
eccentric_orbit_newt.o: eccentric_orbit_newt.c kepler.c mikkola.c
	gcc -Wall -I/home/pzimmerman/include -c eccentric_orbit_newt.c kepler.c mikkola.c
pn_ecc_orbit.o: pn_ecc_orbit.c 
	gcc -Wall -I/home/pzimmerman/include -c pn_ecc_orbit.c pn_kepler.c conservative_dynamics.c
relative_separation.o: relative_separation.c kepler.c
	gcc -Wall -I/home/pzimmerman/include -c relative_separation.c kepler.c
relative_separation_guide.o: relative_separation_guide.c kepler.c
	gcc -Wall -I/home/pzimmerman/include -c relative_separation_guide.c kepler.c
relative_separation_mikkola.o: relative_separation_mikkola.c mikkola.c
	gcc -Wall -I/home/pzimmerman/include -c relative_separation_mikkola.c mikkola.c

clean: 
	rm kepler.o mikkola.o relative_separation_mikkola relative_separation_mikkola.o  relative_separation_guide relative_separation_guide.o 
