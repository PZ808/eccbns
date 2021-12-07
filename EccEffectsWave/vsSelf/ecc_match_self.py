#!/usr/bin/env python

import qm
from optparse import OptionParser
from pylab import sin,cos,sqrt,arange

parser = OptionParser()
parser.add_option("-c", "--con-pnorder", type="int", dest="cpn",
  help="conservative PN order")
parser.add_option("-r", "--rad-pnorder", type="int", dest="rpn",
  help="radiative PN order")   
parser.add_option("-m", "--mass-1", type="float", dest="m1")
parser.add_option("-M", "--mass-2", type="float", dest="m2")
parser.add_option("-f", "--f-min", type="float", dest="fmin")
parser.add_option("-d", "--ligo-type", type="string", dest="detector", 
  help='initial or advanced')

if __name__=='__main__':

  (opts, args) = parser.parse_args()
                               
  eVec  = arange(0.0, 0.41, 0.01)
  f_min = opts.fmin  # fmin causes problems becuase numpy has fmin
  mtot  = opts.m1 + opts.m2
  f_max = qm.schwarzschild_isco(mtot)
  dt    = 1.0 / 4096.0  

  if opts.detector == 'initial': 
    N = 1048576  
    p = qm.new_ligo_psd(N, dt)
    t_max = 256.0
  elif opts.detector == 'advanced':
    N = 16777216  # big for advLIGO
    p = qm.new_adv_ligo_psd(N, dt)
    t_max = 4096.0
  else:
    print 'Error in --ligo-type'

  htilde   = qm.cplx_vec_t(N/2+1)
  stilde   = qm.cplx_vec_t(N/2+1)
  fwd_plan = qm.real_fft_plan_t(N, 1, 0)           
  
  hp0, hc0 = qm.x_model_eccbns_waveform(
             opts.cpn, opts.rpn,
             opts.m1, opts.m2, 0.0, f_min, 0.0, 
             1.0e-16, 0.0, 0.0, 0.0, t_max, 4096.0
             )

  qm.forward_real_fft(htilde, hp0, fwd_plan)

  outfile = open("MatchSelf_"+str(opts.m1)+'_'+str(opts.m2)+'_'+\
str(opts.cpn)+'c'+str(opts.rpn)+'rPN.txt', 'a')

  for ecc in eVec:
    # compute the eccentric waveforms 
    hp, hc = qm.x_model_eccbns_waveform(opts.cpn, 
                                         opts.rpn, 
                                         opts.m1, 
                                         opts.m2, 
                                         ecc, f_min, 
                                         0.0, 1.0e-16,
                                         0.0, 0.0, 
                                         0.0, t_max,
                                         4096.0)         

    qm.forward_real_fft(stilde, hp, fwd_plan)

    m = qm.fd_fd_match(None, stilde, htilde, p, 
      f_min, None, None, None)

    print >> outfile, " %.2f %.6f" %(ecc, m.max())
