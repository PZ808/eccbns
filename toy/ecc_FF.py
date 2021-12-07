#!/usr/bin/python

#
# ecc_FF.py
# 
# Optimizes over chirp mass to compute the FF 
# using a golden search routine 
# Designed to take input from a shell script
#
from pylab import * 
from scipy import optimize
import qm
import math
import sys    

def mismatch(m_chirp):
	global s
	global p
	global amp 
	amp   = qm.sp_amplitude_mc(None,m_chirp,eta,f_min,qm.schwarzschild_isco(m1+m2),N,dt,None) 
    	global phase
    	phase = qm.sp_sstpn_phase_mc(None,m_chirp,eta,0,0,f_min,qm.schwarzschild_isco(m1+m2),N,dt,None) 	
	global htilde 
	htilde = amp*phase
	global m
    	m = qm.td_fd_match(None,s,htilde,p,f_min,None,None,None)
	return 1.0-m.max()

	
if __name__=='__main__':

	if len( sys.argv ) != 4:
		print >> sys.stderr, "Usage: %s filename m1 m2" %sys.argv[0]
		sys.exit(1)  

	in_file = sys.argv[1]
	m1  	= float(sys.argv[2])
	m2  	= float(sys.argv[3])
	m_chirp = float(pow(m1*m2, 3./5.)/pow(m1+m2, 1./5.))

	N 		= 1048576
	f_min	= 40.0
	eta 	= 1.0
	s		= qm.real_vec_t(N)
	s.load(in_file)

	dt    = s.dx
	t     = arange(s.n)*s.dx
	p     = qm.new_ligo_psd(N,dt)

	
	lower 	 = m_chirp
	upper	 = 1.2*m_chirp
	output   = optimize.golden(mismatch,brack=[lower,upper],tol=0.0001,full_output=1)
    	optmass  = output[0]
	optmatch = 1 - output[1]
	outfile  = open('ecc_FF.txt', 'a')
	print >> outfile, "%s:\t FF = %.10f" %(in_file,optmatch)   
	
   
	
	
