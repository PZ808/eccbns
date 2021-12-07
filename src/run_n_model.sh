# $Id: run_n_model.sh,v 1.2 2009/03/07 23:44:00 pzimmerman Exp $
#!/bin/bash
# 
# Automates execution of run_n_model
#

#for e_0 in 0.05 0.10 0.15 0.20 0.25;
#do
#  for M in 2.8 
#  do
#		for f_0 in 40.0 
#		do
#	   	echo "creating ${M}_${f_0}_${e_0} 0-CpN 0-RpN file"  
#	    	./run_n_model ${M} 0.25 ${f_0} ${e_0} 0 0 
#		done
#  done
#done

for e_0 in 0.05 0.10 0.15 0.20 0.25;
do
  for M in 2.8
  do
		for f_0 in 40.0 
		do
	    echo "creating ${M}_${f_0}_${e_0} 2-CpN 0-RpN file"  
	    	./run_n_model ${M} 0.25 ${f_0} ${e_0} 2 0 
		done
  done
done
