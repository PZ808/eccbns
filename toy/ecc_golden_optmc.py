#!/usr/bin/python

#
# ecc_golden_optmc.py
# 
# Optimizes over chirp mass to compute the FF 
# using a golden search routine 
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

	in_file = sys.stdin.next()
	signal  = in_file
	signal 	= signal[0:len(signal)-1]
	m1  	= float(sys.stdin.next())
	m2  	= float(sys.stdin.next())
	m_chirp = float(pow(m1*m2, 3./5.)/pow(m1+m2, 1./5.))

	N 		= 1048576
	f_min	= 40.0
	eta 	= 1.0

	s		= qm.real_vec_t(N)
	s.load(signal) 	

	dt    = s.dx
	t     = arange(s.n)*s.dx
	p     = qm.new_ligo_psd(N,dt)

	
	lower 	 = m_chirp
	upper	 = 1.2*m_chirp
	output   = optimize.golden(mismatch,brack=[lower,upper],tol=0.0001,full_output=1)
    	optmass  = output[0]
	opt_ratio = optmass/m_chirp
	optmatch = 1 - output[1]
	outfile  = open('gold_match.txt', 'a')
	print >> outfile,"%s:\t opt_mc = %.7f, opt_ratio = %.7f, FF = %.7f" %(signal,optmass,opt_ratio,optmatch)   
	
   
	
	
