#!/bin/bash
# $Id: to_dat.sh,v 1.4 2009/06/08 18:07:10 pzimmerman Exp $

# set the sampling interval to 1/4096.0
dt=0.000244140625

eccmodel=xwave

for file in ${eccmodel}*.txt
do 	
  echo "converting ${file}"
  output=`basename $file .txt`.dat
  echo  "% dx = ${dt}" > ${output}
  awk '{print $5}' $file >> ${output}
done
