#
#  run_merc.sh   
#   
#

max_time=$1
max_time=${max_time:?"missing."}
samp=$2
samp=${samp:?"missing."}
max_out=$3
max_out=${max_out:?"missing."}

for t_max in $max_time
do
	for sampling_interval in $samp
  do
  	for out_interval in $max_out
		do
			nohup time ./eccentric_orbit_newt 0.0 ${t_max} 0.0 0.000001 0.000001 ${sampling_interval} ${out_interval} 1> merc.out 3> merc.err < /dev/null &
		done
	done
done

