/********************************************************
 *                                                      *
 * orbit.c                                              *
 *                                                      *
 * Evolves the orbital elements of an eccentric         *
 *        orbit according to radiation reaction         *
 * The orbital variables are eccentricity e,            * 
 *        semilatus rectum p, and orbital phase phi     * 
 *                                                      *
 ********************************************************/
#include <stdio.h>
#include <math.h>
#include <nr.h>
#include <nrutil.h>

#define NUMORBIT 5000   /* number of tabulated points (orbit)  */

const char filename[] = "orbit.txt";  /* data file */

const float mass_1 = 1.4;    /* first mass (solar mass)   */
const float mass_2 = 1.4;    /* second mass (solar mass)  */

const float e_init = 0.2;         /* initial eccentricity */
const float f_init = 30.;         /* initial GW frequency */
const float t_final = 46.14;      /* final time           */
const float t_init = 0.;          /* initial time         */ 

float p[NUMORBIT+1],         /* semilatus (tabulated)     */
      e[NUMORBIT+1],         /* eccentricity (tabulated)  */
      phase[NUMORBIT+1],     /* phase                     */
      f_final,               /* final GW frequency        */
      reduced_mass,          /* reduced mass (in seconds) */
      total_mass;            /* total mass (in seconds)   */

void calc_masses();
void evolve();
void print_info();
void file_results();

int main()
{
  calc_masses();
  evolve();
  print_info();
  file_results();
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

  total_mass = (mass_1 + mass_2)*solar_mass;
  reduced_mass = (mass_1*mass_2/(mass_1+mass_2))*solar_mass;
}

/********************************************************
*                                                       *
* function: print_info                                  *
* Prints relevant quantities to the screen              *
*                                                       *
*********************************************************/
void print_info()
{
  const float PIE = 3.14159265358979;

  f_final = pow( (1.-pow(e[NUMORBIT],2))/p[NUMORBIT], 3./2.)/
            (PIE*total_mass);

  printf("\n");
  printf("mass 1       = %f solar masses\n", mass_1);
  printf("mass 2       = %f solar masses\n", mass_2);
  printf("\n");
  printf("duration     = %f seconds\n", t_final);
  printf("\n");
  printf("initial e    = %f\n", e[1]);
  printf("initial f    = %f Hz\n", f_init);
  printf("initial p    = %f\n", p[1]);
  printf("\n");
  printf("final e      = %f\n", e[NUMORBIT]);
  printf("final f      = %f Hz\n", f_final);
  printf("final p      = %f\n", p[NUMORBIT]);
  printf("\n");
  printf("orbit cycles = %f\n", phase[NUMORBIT]/(2.*PIE));
  printf("orbit points = %d\n", NUMORBIT);
  printf("\n");
}

/********************************************************
*                                                       *
* function: diffeq_orbit                                *
* Defines the differential equations for orbit          *
*                                                       *
*********************************************************/
void diffeq_orbit(float t, float y[], float dydx[])
{
  float p, e, e_squared, mass_factor, e_factor;

  p = y[1];
  e = y[2];
  e_squared = pow(e,2);
  e_factor = pow(1.-e_squared,3./2.);
  mass_factor = (1./reduced_mass)*pow(reduced_mass/total_mass,2);

  dydx[1] = -(64./5.)*mass_factor*(e_factor/pow(p,3))
            *(1.+7.*e_squared/8.);
  dydx[2] = -(304./15.)*mass_factor*e*(e_factor/pow(p,4))
            *(1.+121.*e_squared/304.);
  dydx[3] = pow(1.+e*cos(y[3]),2)/(total_mass*pow(p,3./2.));
}

/********************************************************
*                                                       *
* function: evolve                                      *
*                                                       *
* Prepares the initial values, then integrates          *
* Integration is broken into NUMPHASE intervals         *
* Results from the previous step are used as initial    *
*    data for the next step                             *
* Data is stored in y[1..NVAR][0..NUMPHASE]             *
*                                                       *
*********************************************************/
void evolve()
{
  void init_values(float y_init[]);
  void ode_integrate(void (*derivs)(float, float [], float []),
        float ystart[], int nvar, float x1, 
        float x2, float eps, float h1, float hmin);

  const float PIE = 3.14159265358979;
  const float eps = 1.e-6;
  const float HGUESS = 1.e-3;
  const float HMIN = 0.;

  int count;     /* counts data points               */
  float t1,      /* initial position in current step */
        t2,      /* final position in current step   */
        Delta,   /* length of each interval          */
        *y_init; /* initial values                   */

  Delta = (t_final - t_init)/(NUMORBIT-1);
  y_init = vector(1,3);

  /* preparation of initial values */
  y_init[1] = (1.-e_init*e_init)/pow(PIE*total_mass*f_init,2./3.);
  y_init[2] = e_init;
  y_init[3] = 0.;
  p[1] = y_init[1];
  e[1] = y_init[2];
  phase[1] = y_init[3];

  /* integration */
  for (count = 2; count <= NUMORBIT; count++)
    {
      t1 = t_init + (count-2)*Delta;
      t2 = t1 + Delta;
      ode_integrate(diffeq_orbit,y_init, 3, t1, t2, eps, HGUESS, HMIN);
      p[count] = y_init[1];
      e[count] = y_init[2];
      phase[count] = y_init[3];
    }
}

/********************************************************
*                                                       *
* function: file_results                                *
* Writes the results in a data file                     *
*                                                       *
*********************************************************/
void file_results()
{
  FILE *output_file; 
 
  int count;
  float Delta,t;

  output_file = fopen(filename, "w");

  Delta = (t_final - t_init)/(NUMORBIT-1);

  for (count = 1; count <= NUMORBIT; count++)
    { 
      t = t_init + (count-1)*Delta;
      fprintf(output_file,"%f   %f   %f   %f \n", 
      t,p[count],e[count],phase[count]);
    }

  fclose(output_file);
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

void ode_integrate(void (*derivs)(float, float [], float []), 
     float ystart[], int nvar, float x1, 
     float x2, float eps, float h1, float hmin)
{
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
		(*derivs)(x,y,dydx);
		for (i=1;i<=nvar;i++)
			yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
		if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
		bsstep(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
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



