/* 
 *  reactive_dynamics.c 
 */

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>

double x_dot_0pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( (2.* (37.*e*e*e*e + 292.*e*e + 96.) * eta) /
		(15.* pow(e_factor, 7./2.)) );
}
double e_dot_0pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( -(e*eta * (121.e*e + 304.*e)) /
		(15.* pow(e_factor, 5./2.)) );
}
double x_dot_1pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( (eta / (420.* pow(e_factor, 9./2.))) *
		( -(8288.* eta - 11717.) * pow(e, 6.)
			- 14.* (10122.* eta - 12217.) * pow(e, 4.)
			- 120.* (1330.* eta - 731.) * e*e 
			- 16.* (924.* eta + 743.) ) );
}
double e_dot_1pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( ((e*eta) / (2520.* pow(e_factor, 7./2.))) * 
		( (93184.* eta - 125361.)*e*e*e*e 
			+ 12.* (54271.*eta - 59834.) * e*e 
			+ 8.* (2588.*eta + 8451.) ) );
}
double x_dot_2pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( ( eta / (45360.* pow(e_factor, 11./2.)) ) * 
		( (1964256.* eta*eta - 3259980.* eta + 3523113.) * pow(e, 5.)
			(64828848.* eta*eta - 123108426.* eta + 83424402.) * pow(e, 6.)
			(16650606060.* eta*eta - 20720426.* eta + 783768.) * e*e*e*e
			(61282032.* eta*eta + 15464736.* eta - 92846560.) *e*e +
			1909104.* eta + sqrt(e_factor) * 
			( (264600. - 1058400.* eta) * pow(e, 6.) +
				(64532160. - 25812864.* eta) * e*e - 580608.*eta + 1451520.) +
			4514976.* eta - 360224. ) );
}
double e_dot_2pn ( double e, double eta )
{
	double e_factor = 1. - e*e;
	return ( - ( (e*eta) / (30240.* pow(e_factor, 9./2.)) ) *
		( (2758560.* eta*eta - 4344852.* eta + 3786543.) * pow(e, 6.) +
			(42810096.* eta*eta - 78112266.* eta + 46579718.) * e*e*e*e +
			(48711348.* eta*eta - 35583228.* eta - 36993396.) * e*e +
			4548096.* eta*eta + sqrt(e_factor) *
			((2847600. - 1139040.* eta) * e*e*e*e +
			 (35093520. - 14037408.* eta) *e*e - 5386752.*eta  + 13466880.) +
			13509360. *eta - 15198032. ) );
}

int radiative_evolution( double t, const double y[],
		double f[], void* params )
{
	double dxdt, dedt;

	struct rad_params *rp = (struct rad_params *)params;
	int pn_order = rp->pn_order;
	double eta   = rp->eta;
	double x = y[0];
	double e = y[1];

	double x_fourth  = x*x*x*x;
	double x_fifth   = x*x*x*x*x;
	double x_sixth   = x*x*x*x*x*x;
	double x_seventh = x*x*x*x*x*x*x;

	if ( pn_order == 0 )
	{
		dxdt = x_dot_0pn(e,eta)*x_fifth;
		dedt = e_dot_0pn(e,eta)*x_fourth;
	}
	else if ( pn_order == 1 )
	{
		dxdt = x_dot_0pn(e,eta)*x_fifth + x_dot_1pn(e,eta)*x_sixth;
		dedt = e_dot_0pn(e,eta)*x_fifth + e_dot_1pn(e,eta)*x_sixth;
	}
	else if ( pn_order == 2 )
	{
		dxdt = x_dot_0pn(e,eta)*x_fifth + x_dot_1pn(e,eta)*x_sixth
			+ x_dot_2pn(e,eta)*x_seventh;
		dedt = e_dot_0pn(e,eta)*x_fifth + e_dot_1pn(e,eta)*x_sixth
			+ e_dot_2pn(e,eta)*x_sixth;
	}
	else 
	{
		fprintf( stderr, "Error in pN order: %d\n", pn_order );
		exit(1);
	}

	f[0] = dxdt;
	f[1] = dedt;

	return GSL_SUCCESS;
}

int reactive_dynamics_solver( double abs_step_error
		double rel_step_error, double y_step_scaling,
		double dydt_step_scaling )
{

	const gsl_odeiv_step_type *solver_type
		= gsl_odeiv_step_rkf45;
	gsl_odeiv_step *solver_step
		= gsl_odeiv_step_alloc( solver_type, 2 );
	gsl_odeiv_control *solver_control
		= gsl_odeiv_control_standard_new(
				abs_step_error, rel_step_error,
				y_step_scaling, dydt_step_scaling );
	gsl_odeiv_evolve *solver_evolve
		= gsl_odeiv_evolve_alloc( 2 );
	gsl_odeiv_system solver_system = { 
 		radiative_evolution, NULL, 2
		&ode_params };

	return (gsl_odeiv_evolve_apply(
					solver_evolve, solver_control, 
					solver_step, &solver_system, y) );

}                                    

#if 0
int main ( int argc, char *argv[] )
{
	int i, i_max;
	double t, t_next;
	double step, step_size;
	double dt;

	struct rad_params ode_params;


	dt = 1./4096.;

  ode_params.M 

#endif
