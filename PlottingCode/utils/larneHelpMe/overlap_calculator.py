import sys
import qm
import utils
from numpy import array, arange
from optimize import Amoeba
from pylab import cos

from scipy import interpolate

try:
  import warnings
  warnings.filterwarnings('ignore', module='qm')
except ImportError:
  pass

class overlap_calculator:
  N  = 1048576
  dt = 1.0/4096.0
  
  # allocate the memory needed for match computations
  fwd_plan = qm.real_fft_plan_t(N,1,0)
  rev_plan = qm.cplx_fft_plan_t(N,0,0)
  h        = qm.real_vec_t(N)
  h.dx     = dt
  htilde   = qm.cplx_vec_t(N/2+1)
  stilde   = qm.cplx_vec_t(N/2+1)
  q        = qm.cplx_vec_t(N)
  qtilde   = qm.cplx_vec_t(N)
  overlap  = qm.real_vec_t(N)
  overlap_phase = qm.real_vec_t(N)
  spphase  = qm.cplx_vec_t(N/2+1)
  minus_seven_by_six = qm.new_kfac_vec( N/2+1, -7.0/6.0, 1.0/float(N) )
  minus_one_by_three = qm.new_kfac_vec( N/2+1, -1.0/3.0, 1.0 )

  which_psd = 'initial'
  ligo_psd  = None
  f_low     = -1
  
  phase    = []
  signal_total_mass = 0.0
  pn_approximant    = 'Taylor F2'
  pn_order          = 7

  def __init__(self,waveform_name,m,approximant='Taylor F2',order=7,waveform_eta=0.25,waveform_fc=1024, which_psd='initial'):
    self.signal_total_mass = m
    self.pn_approximant    = approximant
    self.pn_order          = order

    # Load or create the signal
    if waveform_name == 'Taylor 3.5 PN':
      tmp = self.get_htilde(m,waveform_eta,waveform_fc)
      self.stilde.dx = tmp.dx
      for i in range(tmp.n):
        c = qm.cplx_t(tmp[i].re, tmp[i].im)
        self.stilde[i] = c
    elif waveform_name == 'Taylor T4':
      self.pn_approximant = 'Taylor T4'
      self.stilde = self.get_htilde(m, 0.25, qm.schwarzschild_isco(m))
      self.pn_approximant    = approximant
    else:
      # load in the highest resolution as the signal
      s = utils.read_waveform(m/2,waveform_name)

      # For large masses our original N may be too small!
      self.N  = s.n
      self.fwd_plan = qm.real_fft_plan_t(self.N,1,0)
      self.rev_plan = qm.cplx_fft_plan_t(self.N,0,0)
      self.h        = qm.real_vec_t(self.N)
      self.h.dx     = self.dt
      self.htilde   = qm.cplx_vec_t(self.N/2 + 1)
      self.stilde   = qm.cplx_vec_t(self.N/2 + 1)
      self.q        = qm.cplx_vec_t(self.N)
      self.qtilde   = qm.cplx_vec_t(self.N)
      self.overlap  = qm.real_vec_t(self.N)
      self.overlap_phase = qm.real_vec_t(self.N)
      self.spphase  = qm.cplx_vec_t(self.N/2 + 1)
      self.minus_seven_by_six = qm.new_kfac_vec( self.N/2 + 1, -7.0/6.0, 1.0/float(self.N) )
      self.minus_one_by_three = qm.new_kfac_vec( self.N/2 + 1, -1.0/3.0, 1.0 )

      if which_psd == 'initial':
        self.ligo_psd = qm.new_ligo_psd(self.N,self.dt)
        self.f_low    = 40.0

        # self.ligo_psd = qm.new_bench_62_adv_ligo_psd(self.N,self.dt)
        # self.f_low    = 10.0

        # Use the interpolated version
        #self.ligo_psd = self.get_interpolated_psd(self.N,self.dt)
        #self.f_low    = 10.0
        
      else:
        print "Unexpedted PSD requested: ", which_psd
        
      qm.forward_real_fft(self.stilde,s,self.fwd_plan)

  def get_interpolated_psd(self,N,dt):
    bench_f = []
    bench_S = []

    fp = open('/home/lppekows/nr_waveforms/data/AdvLIGO-noise.dat','r')
    for line in fp.xreadlines():
      f,S = line.split('\t')
      bench_f.append(float(f))
      bench_S.append(float(S)*float(S)*qm.DYN_RANGE_FAC*qm.DYN_RANGE_FAC)

    Sn_f = interpolate.interp1d(bench_f,bench_S,bounds_error=False,fill_value=bench_S[0])

    fp.close()

    df = 1.0/(N*dt)
    f = arange(N/2+1) * df
    p = qm.real_vec_t(Sn_f(f))
    p.dx = df

    return p
  



    
  def get_stilde(self):
    return self.stilde

  def get_htilde(self,mass,eta,f_high):
    if self.pn_approximant == 'Taylor F2':
      # Generate a TaylorF2 3.5 frequency domain template
      # which is the same as SPAc(3.5) in Pan et. al.
      qm.sp_amplitude_mc(self.htilde,      # array where result is stored
                         mass,             # Mass parameter
                         eta,              # eta parameter
                         self.f_low,       # starting frequency
                         f_high,           # ending frequency
                         self.N,           # number of points to calc
                         self.dt,   
                         self.minus_seven_by_six)
      
      qm.sp_sstpn_phase_mc(self.spphase,   # Result array
                           mass,
                           eta, 
                           0,              # beta (used for 1.5 pn)
                           self.pn_order,  # order * 2 (eg, 7/2 = 3.5)
                           self.f_low,     # start freq
                           f_high,         # end freq
                           self.N,
                           self.dt,
                           self.minus_one_by_three)
      
      qm.mul_cplx_vec_in_place(self.htilde,self.spphase)
    elif self.pn_approximant == 'Taylor T4':
      a_plus, a_cross, phi = qm.tpn_waveform(mass, eta, self.f_low,
                                             self.N, self.dt)

      # Phi    = 2*phi
      # h_stpn = a_plus.array()*cos(Phi)

      for i in range(len(a_plus)):
        a_plus[i] = a_plus[i] * cos(2 * phi[i])

      # Go to freq domain
      qm.forward_real_fft(self.htilde, a_plus, self.fwd_plan)
      
    return self.htilde


  def calc_overlap(self,mass,eta,f_high):
    ht = self.get_htilde(mass,eta,f_high)

    # Now calculate the match
    qm.fd_fd_match(self.overlap,
                   self.stilde,
                   ht,
                   self.ligo_psd,
                   self.f_low,
                   self.q,
                   self.qtilde,
                   self.rev_plan)
      
    return self.overlap, self.phase


  def calc_overlap_max(self,mass,eta,f_high):
    if (mass < 0.0) or (f_high > 2048):
      return -1

    overlap, phase = self.calc_overlap(mass,eta,f_high)
    tmp = overlap.max()

    # This is a kludge: sometimes this returns 'inf.'  That only seems
    # to happen for masses that are well-outside the band, but...

    if tmp > 2.0:
      print "Death: ", mass, eta, f_high
      return -1
      
    return tmp


  def mass_freq_grid(self,m_low=None,m_high=None,m_steps=15,
                     f_low=100,f_high=8196,f_steps=15,eta=0.25):

    if m_low == None:
      m_low = self.signal_total_mass * (0.5)

    if m_high == None:
      m_high = self.signal_total_mass * (1.5)

    m_inc = (m_high-m_low)/m_steps
    xs    = map((lambda x: x * m_inc + m_low),range(0,m_steps+1))

    f_inc = (f_high-f_low)/f_steps
    ys    = map((lambda x: x * f_inc + f_low),range(0,f_steps+1))

    x_grid = []
    y_grid = []
    z_grid = []

    for x in xs:
      x_row = []
      y_row = []
      z_row = []

      for y in ys:
        x_row.append(x)
        y_row.append(y)
        z_row.append(self.calc_overlap_max(x,eta,y))
        
      x_grid.append(x_row)
      y_grid.append(y_row)
      z_grid.append(z_row)

    return x_grid,y_grid,z_grid


  def optimize_over_template_bank(self,f_c_min = 100, f_c_max = 2048, f_steps = 25):
    best_overlap = previous_best_overlap = -1
    best_m   = 0
    best_f_c = 0
    best_eta = 0

    m_low  = self.signal_total_mass * 0.5
    m_high = self.signal_total_mass * 2

    f_min   = 100
    f_max   = 1500
    f_steps = 15

    for i in [0,1,2]:
      f_inc    = (f_max-f_min)/f_steps
      f_highs  = map((lambda x: x * f_inc + f_min),range(f_steps+1))
      previous = -1
      
      for f_high in f_highs:
        if f_high == previous or f_high < 50:
          continue

        sys.stdout.flush()
        
        previous = f_high
        
        n, m1, m2, t0, t3 = qm.new_lal_bank( m_low,         # Min mass
                                             m_high,        # Max mass
                                             0.97,          # Min match
                                             self.ligo_psd, # PSD
                                             self.f_low,    # f_c min
                                             f_high)        # f_c max

        for i in range(n):
          m_total = m1[i] + m2[i]
          eta     = m1[i] * m2[i] / (m_total * m_total)

          overlap = self.calc_overlap_max(m_total,eta,f_high)

          if overlap > best_overlap:
            best_overlap = overlap
            best_m       = m_total
            best_eta     = eta
            best_f_c     = f_high

        f_min = best_f_c - f_inc
        f_max = best_f_c + f_inc
          
    return best_m, best_eta, best_f_c, best_overlap

  def optimize_over_f_c(self, mass, eta, f_min = 100, f_max = 2048, f_steps = 15):
    best_overlap = previous_best_overlap = -1
    best_f_c = 0


    for i in [0,1,2]:
      f_inc    = (f_max-f_min)/f_steps
      f_highs  = map((lambda x: x * f_inc + f_min),range(f_steps+1))
      previous = -1
      
      for f_high in f_highs:
        if f_high == previous or f_high < 50:
          continue

        previous = f_high        
        overlap  = self.calc_overlap_max(mass,eta,f_high)

        if overlap > best_overlap:
          best_overlap = overlap
          best_f_c     = f_high

      f_min = best_f_c - f_inc
      f_max = best_f_c + f_inc
          
    return best_f_c, best_overlap
    

  def find_error_bars_helper(self, f, start_value, start_overlap, drop, incr, low, high):
    upper_range   = []
    upper_indices = []
    lower_range   = []
    lower_indices = []
    
    current = start_value
    overlap = start_overlap
      
    while start_overlap - overlap < drop and current < high:
      current += incr
      overlap  = f(current)
      upper_indices.append(current)
      upper_range.append(overlap)

    current = start_value
    overlap = start_overlap

    while start_overlap - overlap < drop and current > low:
      current -= incr
      overlap  = f(current)
      lower_indices.append(current)
      lower_range.append(overlap)

    lower_range.reverse()
    lower_indices.reverse()
      
    full_range   = lower_range + [start_overlap] + upper_range
    full_indices = lower_indices + [start_value] + upper_indices
    
    return full_indices, full_range

  
  def find_m_error_bars(self, m, eta, f_c, drop, low, high):
    # Find the starting overlap
    start_overlap = self.calc_overlap_max(m, eta, f_c)
    f             = lambda m: self.calc_overlap_max(m, eta, f_c)
    return self.find_error_bars_helper(f, m, start_overlap, drop, 0.001, low, high)

  def find_eta_error_bars(self, m, eta, f_c, drop, low, high):
    start_overlap = self.calc_overlap_max(m, eta, f_c)
    f             = lambda eta: self.calc_overlap_max(m, eta, f_c)
    return self.find_error_bars_helper(f, eta, start_overlap, drop, 0.001, low, high)

  def find_f_c_error_bars(self, m, eta, f_c, drop, low, high):
    start_overlap = self.calc_overlap_max(m, eta, f_c)
    f             = lambda f_c: self.calc_overlap_max(m, eta, f_c)
    return self.find_error_bars_helper(f, f_c, start_overlap, drop, 1.0, low, high)


    
  def optimize(self):
    f    = lambda x: 1.0-self.calc_overlap_max(x[0],x[1],x[2])
    isco = qm.schwarzschild_isco(self.signal_total_mass)

    am   = Amoeba(3,
                  f,
                  points,
                  tol=1.0e-4, 
                  lower_limits = [self.signal_mass * 0.5, 0.0, 100],
                  upper_limits = [self.signal_mass * 2.0, 1.0, 1024])

    am.optimize()
    am.restart()
    am.optimize()
    
    
    return res

