# evolve_eccbns.sh   
# $Id: evolve_eccbns.sh,v 1.8 2009/06/08 18:07:10 pzimmerman Exp $  

# use the x_drive driver to run
# simulate the evolution of 
# of an eccentric binary to 
# 3pN in the conservative dynamics (equations of motion)
# and 2pN in the reactive dynamics (GW flux) 

driver=x_drive

# parse arguments
#cPN=$1      # conservative pN order
#cPN=${cPN:?"missing."}
rPN=$1      # reactive pN order
rPN=${rPN:?"missing."}
f_init=$2   # initial GW frequency (f_init <= 40.0 Hz)
f_init=${f_init:?"missing."}
l_init=$3   # initial mean-anomaly (default should be zero)
l_init=${l_init:?"missing."}
ode_eps=$4  # tolerance of the GSL ode driver (default is 1e-16)
ode_eps=${ode_eps:?"missing."}
samp=$5    # samping rate (default is 1./4096.)
samp=${samp:?"missing."}

# Loop over initial eccentricity in 
# increments of 0.05 
for cPN in 1 2 3 
do
  for e_init in $(seq 0.00 0.01 0.1)  $(seq 0.12 0.02 0.40);
  do
    # Loop over masses
    for m1 in 1.4; 
    do
      for m2 in 1.4 5.0 10.0;  
      do
        echo "${driver} evolving ${m1} ${m2} ${e_init} system at ${cPN} con-pN ${rPN}/2 rad-pN"   
        ./${driver} ${cPN} ${rPN} ${m1} ${m2} ${f_init} ${e_init} ${l_init} ${ode_eps} ${samp}
      done
    done
  done
done
