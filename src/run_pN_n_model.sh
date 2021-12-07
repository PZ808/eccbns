# $Id: run_pN_n_model.sh,v 1.6 2009/03/07 23:44:00 pzimmerman Exp $
#!/bin/bash
# 
# Automates execution of evolve2Rad2ConPN
#

for e_0 in 0.05 0.10 0.15 0.20 0.25;
do
  for M in 2.8 
  do
		for f_0 in 40.0 
		do
	   	echo "creating ${M}_${f_0}_${e_0} 0-CpN 0-RpN file"  
	    	./n_model_evolution ${M} 0.25 ${f_0} ${e_0} 0 0 0.000001
		done
  done
done

for e_0 in 0.05 0.10 0.15 0.20 0.25;
do
  for M in 2.8
  do
		for f_0 in 40.0 
		do
	    echo "creating ${M}_${f_0}_${e_0} 2-CpN 0-RpN file"  
	    	./n_model_evolution ${M} 0.25 ${f_0} ${e_0} 2 0 0.000001
		done
  done
done

#nohup time ./n_model_evolution ${M} 0.25 ${f_0} ${e_0} 2 0 1> n_model_2c0r_waveforms.out 2> n_model_2c0r_waveforms.err < /dev/null &
#nohup time ./n_model_evolution ${M} 0.25 ${f_0} ${e_0} 0 0 1> n_model_0c0r_waveforms.out 2> n_model_0c0r_waveforms.err < /dev/null &
