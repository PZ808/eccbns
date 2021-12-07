#!/bin/bash

# compute the match between a T4
# PN order templateOrder against an x-model eccentric 
# signal at orders radOrder,conOrder using compute_match.py

fLow=40.0

# Get the command line arguments 
# pn order of the template
templateOrder=$1
templateOrder=${templateOrder:?"PN order of the template is missing."}
# conservative pn order of the siganl
conOrder=$2
conOrder=${conOrder:?"conPN order of the signal is missing."}
# radiation reaction pn order of the siganl
radOrder=$3
radOrder=${radOrder:?"radPN order of the signal is missing."}
# noise PSD of the detector
#noise=$4
#noise=${noise:?"noise PSD is missing. Use initial for iLIGO or advanced for advLIGO."}

# get the lower frequency cut off
# if [ ${noise} = initial ]; then
#    fLow=40.0
# elif [ ${noise} = advanced ]; then
#     fLow=10.0
# else 
#     echo "$noise is not a valid argument for LIGO detector type"
#     exit 1
# fi

#case ${noise} in
# initial LIGO
#  initial ) 
#    fLow=40.0 ;;
# advanced LIGO
#  advanced )
#    fLow=10.0 ;;
# unknown
#  * ) echo "$noise is not a valid argument for LIGO detector type." 
#      exit 1 ;;
#esac

# loop over the masses we want and compute the overlaps
for m1 in 1.4
do
  for m2 in 1.4 5.0 10.0
  do
    python t4_vs_xmod_match.py -o $templateOrder -c $conOrder -r $radOrder\
 -m $m1 -M $m2 -f $fLow  
  done
done


for m1 in 5.0
do
  for m2 in 5.0
  do
    python t4_vs_xmod_match.py -o $templateOrder -c $conOrder -r $radOrder\
 -m $m1 -M $m2 -f $fLow  
  done
done




