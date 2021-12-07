/********************************************************
 *                                                      *
 * duration.c                                           *
 *                                                      *
 * Evolves the orbital elements of an eccentric         *
 *        orbit according to radiation reaction         *
 * Calculates the duration of the orbit, from a         *
 *        provided initial mean GW frequency to         *
 *        a final value of p                            *
 *                                                      *
 ********************************************************/
#include <stdio.h>
#include <math.h>
#include <nr.h>
#include <nrutil.h>

#define NUMPOINTS 1000   /* number of data points         */
#define NVAR 2           /* number of dependent variables */
#define EPS 1.e-6        /* accuracy parameter            */
#define HGUESS -1.e-3    /* first guess for step size     */
#define HMIN 0.          /* minimum step size             */

const float mass_1 = 1.4;    /* first mass (solar mass)      */
const float mass_2 = 10.;    /* second mass (solar mass)     */

const float e_init = 0.5;    /* initial eccentricity         */
const float f_init = 26.6666; /* initial GW frequency (in Hz) */ 
const float p_final = 6.;     /* final value of p              */

float y[NVAR+1][NUMPOINTS+1], /* dependent variables       */
      p_init,                 /* initial value of p        */
      f_final,                /* final GW frequency        */
      t_final,                /* final value of time       */
      e_final,                /* final eccentricity        */
      reduced_mass,           /* reduced mass (in seconds) */
      total_mass;             /* total mass (in seconds)   */

void calc_masses();
void evolve();
void unpack();
void print_info();

int main()
{
  calc_masses();
  evolve();
  unpack();
  print_info();
  return (0);
}

/********************************************************
*                                                       *
* function: calc_masses                                 *
* Calculates the reduced and total masses in seconds    *
*                                                       *
*********************************************************/
void calc_masses()
{
  const float solar_mass = 4.9255e-6;  /* solar mass in sec */
  const float PIE = 3.14159265358979;

  total_mass = (mass_1 + mass_2)*solar_mass;
  reduced_mass = (mass_1*mass_2/(mass_1+mass_2))*solar_mass;

  p_init = (1.-e_init*e_init)/pow(PIE*total_mass*f_init,2./3.);
}

/********************************************************
*                                                       *
* function: unpack                                      *
* Calculates orbital parameters from y[ ]               *
*                                                       *
*********************************************************/
void unpack()
{
  const float PIE = 3.14159265358979;

  t_final = y[1][NUMPOINTS];
  e_final = y[2][NUMPOINTS];
  f_final = pow( (1.-pow(e_final,2))/p_final, 3./2.)/
            (PIE*total_mass);
}

/********************************************************
*                                                       *
* function: print_info                                  *
* Prints relevant quantities to the screen              *
*                                                       *
*********************************************************/
void print_info()
{
  printf("\n");
  printf("mass 1      = %f solar masses\n", mass_1);
  printf("mass 2      = %f solar masses\n", mass_2);
  printf("\n");
  printf("initial e   = %f\n", e_init);
  printf("initial f   = %f Hz\n", f_init);
  printf("initial p   = %f\n", p_init);
  printf("final e     = %f\n", e_final);
  printf("final f     = %f Hz\n", f_final);
  printf("final p     = %f\n", p_final);
  printf("\n");
  printf("duration    = %f seconds\n", t_final);
  printf("\n");
}

/********************************************************
*                                                       *
* function: init_values                                 *
* Calculates the initial values                         *
     of the dependent variables                         *
*                                                       *
*********************************************************/
void init_values(float y_init[])
{
  y_init[1] = 0.;
  y_init[2] = e_init;
}

/********************************************************
*                                                       *
* function: diff_eqns                                   *
* Defines the differential equations                    *
*                                                       *
*********************************************************/
void diff_eqns(float p, float y[], float dydx[])
{
  float e_simple, e_squared, e_factor, mass_factor;

  e_simple = y[2];
  e_squared = pow(e_simple,2);
  e_factor = pow(1.-e_squared,-3./2.);
  mass_factor = reduced_mass*pow(total_mass/reduced_mass,2);

  dydx[1] = -(5./64.)*mass_factor*pow(p,3)*e_factor/
             (1.+7.*e_squared/8.);
  dydx[2] = (19./12.)*(e_simple/p)*(1.+121.*e_squared/304.)/
            (1.+7.*e_squared/8.);
}

/********************************************************
*                                                       *
* function: evolve                                      *
*                                                       *
* Prepares the initial values, then integrates          *
* Integration is broken into NUMPOINTS intervals        *
* Results from the previous step are used as initial    *
*    data for the next step                             *
* Data is stored in y[1..NVAR][0..NUMPOINTS]            *
*                                                       *
*********************************************************/
void evolve()
{
  void init_values(float y_init[]);
  void ode_integrate(float ystart[], int nvar, float x1, 
        float x2, float eps, float h1, float hmin);
  void diff_eqns(float x, float y[], float dydx[]);

  int n,         /* counts variables                 */
      count;     /* counts data points               */
  float p1,      /* initial position in current step */
        p2,      /* final position in current step   */
        Delta,   /* length of each interval          */
        *y_init; /* initial values                   */

  Delta = (p_final - p_init)/(NUMPOINTS-1);
  y_init = vector(1,NVAR);

  /* preparation of initial values */
  init_values(y_init);  
  for (n = 1; n <= NVAR; n++)
  {
     y[n][1] = y_init[n];      
  }

  /* integration */
  for (count = 2; count <= NUMPOINTS; count++)
    {
      p1 = p_init + (count-2)*Delta;
      p2 = p1 + Delta;
      ode_integrate(y_init, NVAR, p1, p2, EPS, HGUESS, HMIN);
      for (n = 1; n <= NVAR; n++)
      {
         y[n][count] = y_init[n];
      }
    }
}

/********************************************************
*                                                       *
* function: ode_integrate                               *
*                                                       *
* Modified version of NR odeint                         *
*                                                       *
*********************************************************/
#define MAXSTP 10000
#define TINY 1.0e-30

void ode_integrate(float ystart[], int nvar, float x1, 
     float x2, float eps, float h1, float hmin)
{
   void bsstep1(float y[], float dydx[], int nv, float *xx, float htry, 
   float eps, float yscal[], float *hdid, float *hnext,
   void (*derivs)(float, float [], float []));

	int nstp,i;
	float x,hnext,hdid,h;
	float *yscal,*y,*dydx;

	yscal=vector(1,nvar);
	y=vector(1,nvar);
	dydx=vector(1,nvar);
	x=x1;
	h=SIGN(h1,x2-x1);

	for (i=1;i<=nvar;i++) y[i]=ystart[i];

	for (nstp=1;nstp<=MAXSTP;nstp++) {
		diff_eqns(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		bsstep(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,diff_eqns);
		if ((x-x2)*(x2-x1) >= 0.0) {
			for (i=1;i<=nvar;i++) ystart[i]=y[i];
			free_vector(dydx,1,nvar);
			free_vector(y,1,nvar);
			free_vector(yscal,1,nvar);
			return;
		}
		if (fabs(hnext) <= hmin) 
                   nrerror("Step size too small in odeint");
		h=hnext;
	}
	nrerror("Too many steps in routine odeint");
}

