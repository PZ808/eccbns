/* mass of the sun in kg */
#define M_SUN_KG 1.9884e30
/* mass of the sun in seconds */
#define M_SUN_S 4.92549095e-6 
/* mass of the sun in meters */
#define M_SUN_M 1.47662504e3
/* Newton's gravitational constant */
#define G_NEWT 6.6743e-11
/* speed of light */
#define C_SI 299792458.0

struct ode_parameters
{
	double eta;
	int conservative_pn_order;
	int radiation_pn_order;
};

struct kepler_params{
  int pn_order;
  double eta;
  double x;
  double e;
  double l;
};

int eccentric_x_model_odes(
    double t, 
    const double y[],
    double dydt[], 
    void *params );

double dx_dt( int radiation_pn_order, double eta, double x, double e );
double de_dt( int radiation_pn_order, double eta, double x, double e );
double dl_dt( int conservative_pn_order, double eta, double x, double e ); 
double dphi_dt( int conservative_pn_order, double u, double eta, double x, double e );

typedef enum 
{
  kappa_e,
  kappa_j
} kappa_type;   

double kappa( kappa_type kappa_t, int p_max, double e );

double x_dot_0pn( double e, double eta );
double x_dot_1pn( double e, double eta );
double x_dot_1p5pn( double e, double eta );
double x_dot_2pn( double e, double eta );

double e_dot_0pn( double e, double eta );
double e_dot_1pn( double e, double eta );
double e_dot_1p5pn( double e, double eta );
double e_dot_2pn( double e, double eta );

double l_dot_1pn( double e, double eta );
double l_dot_2pn( double e, double eta );
double l_dot_3pn( double e, double eta );

double phi_dot_0pn( double e, double eta, double u );
double phi_dot_1pn( double e, double eta, double u );
double phi_dot_2pn( double e, double eta, double u );
double phi_dot_3pn( double e, double eta, double u );

double pn_kepler_equation( int conservative_pn_order, double eta, double x, double e, double l );
double kepler_equation( double u, void* params );
double mikkola_finder( double e, double l);

double separation( int conservative_pn_order, double u, double eta, double x, double e);

double x_cartesian( double r, double phi );
double y_cartesian( double r, double phi );
