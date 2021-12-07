#!/usr/bin/env python
# $Id: t4_vs_xmod_match.py,v 1.2 2009/07/17 16:58:03 pzimmerman Exp $

#
# t4_vs_xmod_match.py
#
# Script computes the overlap
# between a taylor T4 waveform and 
# and an eccentric x-model waveform.
#

import sys
import qm
from pylab import arange, sin, cos, sqrt
from optparse import OptionParser

# Parse the command line arguments
parser = OptionParser()
parser.add_option("-o", "--t4-pnorder", type="int", dest="tpn",
    help="Taylor T4 waveform PN order")
parser.add_option("-c", "--con-pnorder", type="int", dest="cpn",
    help="Eccentric conservative PN order")
parser.add_option("-r", "--rad-pnorder", type="int", dest="rpn",
    help="Eccentric waveform radiative PN order")   
parser.add_option("-m", "--mass-1", type="float", dest="mass1")
parser.add_option("-M", "--mass-2", type="float", dest="mass2")
parser.add_option("-f", "--f-min", type="float", dest="fLow")

(opts, args) = parser.parse_args()

tpn   = opts.tpn
cpn   = opts.cpn
rpn   = opts.rpn
m1    = opts.mass1
m2    = opts.mass2
f_low = opts.fLow

# Common parameters
###################
fs    = 4096.0   # sampling rate (Hz)
dt    = 1.0/fs   # sampling interval (sec)
N     = 1048576  # number of sample points
iota  = 0.0      # inclination angle
beta  = 0.0      # inclination angle
ell   = 0.0      # mean anomly
toler = 1.0e-16  # tolerance of the ode solver for ecc code 
t_min = 0.0  
t_max = 256.0    # t_max/dt give the value of h.dx for the ecc waveform
m_tot = m1+m2    # total mass
eta   = m1*m2/(m_tot**2)  # symmetic mass ratio

# create a file for writing out the match data
outfile = open('Match_T4vsXmodel_' + str(tpn) + '_' +  
              str(m1) + '_' + str(m2) + '_' + 
              str(cpn) + 'c' + str(rpn) + 'r' + 'PN.txt', 'a') 

# list of eccentricities
ecc = arange(0.00, 0.11, 0.01)
# list for match results
res = []

ligo_psd = qm.new_ligo_psd(N,dt)    
hp_tpn   = qm.real_vec_t(N)
stilde   = qm.cplx_vec_t(N/2 + 1)     
fwd_plan = qm.real_fft_plan_t(N,1,0)

overlap  = qm.real_vec_t(N)         

# Generate the tpn waveform which is a T4 approximant
# in the time domain
a_plus, a_cross, phi = qm.tpn_waveform(m_tot, eta, tpn, f_low, N, dt)                      

t_T4 = arange(0, a_plus.n) * a_plus.dx 
# get the dominant GW harmonic 
Phi  = 2*phi
# put the GW signal for T4 in a numpy array
h_T4 = a_plus.array()*cos(Phi)  

# populate the qm vector with the elements of the 
# numpy array  
hp_tpn.dx = dt
for i in range(hp_tpn.n):
  hp_tpn[i] = a_plus[i]*cos(Phi[i])

#
# Compute the overlap between several 
# eccentric wavforms and the T4 waveform.
#
for e in ecc:
  # compute the eccentric waveforms
  hp, hc = qm.x_model_eccbns_waveform(cpn, rpn, m1, m2, 
    e, f_low, ell, toler, iota, beta, t_min, t_max, fs) 

  # FFT the eccentric waveform into the frequency domain
  qm.forward_real_fft(stilde, hp, fwd_plan)
  
  # compute the match between the T4 and the x-model waveforms for a
  # given eccentricity 
  overlap = qm.td_fd_match(overlap, 
                           hp_tpn, 
                           stilde, 
                           ligo_psd, 
                           f_low, 
                           None, None, None)

  # put results in a list to print to file
  res.append((e, overlap.max()))

# write results to file
for line in res:
  print >> outfile, "%.2f %.4f" %(line)
