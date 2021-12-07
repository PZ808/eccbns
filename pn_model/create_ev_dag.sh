# $Id: create_ev_dag.sh,v 1.4 2009/06/08 18:07:10 pzimmerman Exp $
# create_ev_dag.sh                      
# Generates a condor ready DAG file to run the 
# the x-model code using a cluster             

driver=x_drive  

# parse arguments
cPN=$1      # conservative pN order
cPN=${cPN:?"missing."}
rPN=$2      # reactive pN order
rPN=${rPN:?"missing."}
f_init=$3   # initial GW frequency (f_init <= 40.0 Hz)
f_init=${f_init:?"missing."}
l_init=$4   # initial mean-anomaly (default should be zero)
l_init=${l_init:?"missing."}
ode_eps=$5  # tolerance of the GSL ode driver (default is 1e-16)
ode_eps=${ode_eps:?"missing."}
samp=$6    # samping rate (default is 1./4096.)
samp=${samp:?"missing."}

count=0
Mmax=15.0

for e_init in $(seq 0.00 0.1 0.4);  
do
  for m1 in $(seq 1.0 0.5 15.0); 
  do
    for m2 in $(seq 1.0 0.5 $m1);
    do
      Mtot=`echo $m1+$m2 | bc` # compute total mass as float 
      mtest=`echo "${Mtot} ${Mmax}" | awk '{if ($1 <= $2) print 1; else print 0}'`
      if [ $mtest -ne 0 ]; then
        echo "JOB xdrive${count} xdrive.sub" >> ecc_evolve_${cPN}c${rPN}rPN_005.dag
        echo "VARS xdrive${count} c=\"${cPN}\" r=\"${rPN}\" m=\"${m1}\" q=\"${m2}\" f=\"${f_init}\" e=\"${e_init}\" l=\"${l_init}\" eps=\"${ode_eps}\" dt=\"${samp}\"" >> ecc_evolve_${cPN}c${rPN}rPN_005.dag
        count=$((count+1))
      fi
    done
  done
done
# Alternate way of testing using float
#d=`echo $Mmax - $Mtot |bc`
#val=`echo $d |grep "-" |wc -l`
#if [ $val -ge 0 ]; then
