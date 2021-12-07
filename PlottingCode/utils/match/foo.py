class GoldenFF:
    N=pow(2,20)
    fs=4096.0
    dt=1./fs
    fwd_plan = qm.real_fft_plan_t(N,1,0)
    h        = qm.real_vec_t(N)
    h.dx     = dt
    htilde   = qm.cplx_vec_t(N/2+1)
    stilde   = qm.cplx_vec_t(N/2+1)
    def __init__(self, con_order,rad_order,m1,m2,e0,f0,t_max, psd='initial'):
        self.ecc_signal_con_order = con_order
        self.ecc_signal_rad_order = rad_order
        self.ecc_signal_m1 = m1
        self.ecc_signal_m2 = m2
        tmp_ecc_hp, tmp_ecc_hc = self.get_ecc_signal(con_order,rad_order,m1,m2,e0,f0,t_max)
    def get_ecc_signal(self,con_order,rad_order,m1,m2,ecc,freq_low,ell=0,eps=1.e-16,iota=0,beta=0,t0=0,tmax=256.,fs=4096.):
        return qm.x_model_eccbns_waveform(con_order,rad_order,m1,m2,ecc,freq_low,ell,eps,iota,beta,t0,tmax,fs)

if __name__=='__main__':

  x = GoldenFF(3,4,1.4,1.4,0.0,40.0,256.0)
  t = arange(x.N)*x.dt
  hp, hc = x.get_ecc_signal(3,4,1.4,1.4,0.1,40.0)
  plot(t,hp)
