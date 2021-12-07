struct t4_params
{
  double mass1;
  double mass2;
};

int t4_ode_system( 
    double t, 
    const double y[], 
    double dydt[], 
    void *params 
    );

int tpn_waveform( real_vec_t** h_plus, 
                  real_vec_t** h_cross, 
                  double m1,
                  double m2,
                  double f_min, 
                  int N, 
                  double dt 
                  );            
