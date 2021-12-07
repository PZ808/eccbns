/********************************************************
 *                                                      *
 * wave.c                                               *
 *                                                      *
 * Calculates the Newtonian waveform of an evolving     *
 *        eccentric orbit                               * 
 * Orbital evolution is handled by embedded code        *
 *        from orbit.c                                  *
 *                                                      *
 ********************************************************/
#include <stdio.h>
#include <math.h>
#include <nr.h>
#include <nrutil.h>

#define NUMPOINTS 50000 /* number of data points (wave)  */
#if 1
#undef M_PI
#define M_PI 3.14159265358979
#endif

const char filename[] = "wave.txt";  /* data file */

const float mass_1 = 1.4;    /* first mass (solar mass)   */
const float mass_2 = 1.4;    /* second mass (solar mass)  */

const float e_init = 0.5;         /* initial eccentricity */
const float t_final = 26.22;      /* final time           */

const float t_init = 0.;          /* initial time         */ 
const float f_init = 26.6666;     /* initial GW frequency */

const float iota = 45.;      /* source's first polar angle    */
const float bet = 0.;        /* source's second polar angle   */
const float theta = 45.;     /* detector's first polar angle  */
const float phi = 20.;       /* detector's second polar angle */
const float psi = 0.;        /* detector's third polar angle  */

float p[NUMPOINTS+1],        /* semilatus (tabulated)     */
      e[NUMPOINTS+1],        /* eccentricity (tabulated)  */
      phase[NUMPOINTS+1],    /* phase                     */
      h[NUMPOINTS+1],        /* wave                      */
      f_final,               /* final GW frequency        */
      reduced_mass,          /* reduced mass (in seconds) */
      total_mass;            /* total mass (in seconds)   */

void calc_masses();
void evolve();
void calc_wave();
void print_info();
void file_results();

int main()
{
  calc_masses();
  evolve();
  calc_wave();
  print_info();
  file_results();
  return (0);
}

/********************************************************
* function calc_wave                                    *
* Calculates the gravitational waves                    *
*                                                       *
*********************************************************/
void calc_wave()
{
  const float PIE = M_PI;

  int count;
  float pp,ee,pha,normalization,t,Delta;
  float c,c2,s,s2,iotaa,betaa,thetaa,phii,psii;
  float h_plus,h_cross,f_plus,f_cross;

  iotaa = iota*PIE/180.;
  betaa = bet*PIE/180;
  c = cos(iotaa);
  c2 = pow(c,2);
  s = sin(iotaa);
  s2 = pow(s,2);
  thetaa = theta*PIE/180.;
  phii = phi*PIE/180.;
  psii = psi*PIE/180.;
  normalization = p[1];

  f_plus = 0.5*(1.+pow(cos(thetaa),2))*cos(2.*phii)*cos(2.*psii)
           - cos(thetaa)*sin(2.*phii)*sin(2.*psii);
  f_cross = 0.5*(1.+pow(cos(thetaa),2))*cos(2.*phii)*sin(2.*psii)
            + cos(thetaa)*sin(2.*phii)*cos(2.*psii);

  Delta = (t_final - t_init)/(NUMPOINTS-1);
  fprintf(stderr, "delta = %32.16e\n", Delta );

  for (count = 1; count <= NUMPOINTS; count++)
    {
      t = t_init + (count-1)*Delta;
      pp = p[count];
      ee = e[count];
      pha = phase[count];

      h_plus = ( 2.*(1.+c2)*cos(2.*(pha-betaa))
         + (ee/2.)*(1+c2)*(5.*cos(pha-2.*betaa)
                          + cos(3.*pha-2.*betaa))
         + ee*s2*cos(pha)
         + pow(ee,2)*((1.+c2)*cos(2.*betaa)+s2) )
         *normalization/pp;

      h_cross = ( 4.*c*sin(2.*(pha-betaa))
         + ee*c*(5.*sin(pha-2.*betaa)
                + sin(3.*pha-2.*betaa))
         - 2.*pow(ee,2)*c*sin(2.*betaa) )
         *normalization/pp;

      h[count] = f_plus*h_plus + f_cross*h_cross;
    }      
}

/********************************************************
*                                                       *
* function: calc_masses                                 *
* Calculates the reduced and total masses in seconds    *
*                                                       *
*********************************************************/
void calc_masses()
{
#if 1
  /* value used by Poisson and Martel */
  const float solar_mass = 4.9255e-6;  /* solar mass in sec */
#else
  /* value used in LIGO */
  const float solar_mass = 4.92549095e-06;  /* solar mass in sec */
#endif

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
  const float PIE = M_PI;

  f_final = pow( (1.-pow(e[NUMPOINTS],2))/p[NUMPOINTS], 3./2.)/
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
  printf("final e      = %f\n", e[NUMPOINTS]);
  printf("final f      = %f Hz\n", f_final);
  printf("final p      = %f\n", p[NUMPOINTS]);
  printf("\n");
  printf("orbit cycles = %f\n", phase[NUMPOINTS]/(2.*PIE));
  printf("wave points  = %d\n", NUMPOINTS);
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
* Integration is broken into NUMPOINTS intervals        *
* Results from the previous step are used as initial    *
*    data for the next step                             *
* Data is stored in y[1..NVAR][0..NUMPOINTS]            *
*                                                       *
*********************************************************/
void evolve()
{
  void init_values(float y_init[]);
  void ode_integrate(void (*derivs)(float, float [], float []),
        float ystart[], int nvar, float x1, 
        float x2, float eps, float h1, float hmin);

  const float PIE = M_PI;
  const float eps = 1.e-6;
  float HGUESS = 1.e-3;
  const float HMIN = 0.;

  int count;     /* counts data points               */
  float t1,      /* initial position in current step */
        t2,      /* final position in current step   */
        p_init,
        Delta,   /* length of each interval          */
        *y_init; /* initial values                   */

  Delta = (t_final - t_init)/(NUMPOINTS-1);
  fprintf(stderr, "delta = %32.16e\n", Delta );
  y_init = vector(1,3);
  p_init = (1.-e_init*e_init)/pow(PIE*total_mass*f_init,2./3.);
  fprintf( stderr, "PIE        = %32.16e\n", PIE );
  fprintf( stderr, "total_mass = %32.16e\n", total_mass );
  fprintf( stderr, "f_init     = %32.16e\n", f_init );
  fprintf( stderr, "p_init     = %32.16e\n", p_init );
  fprintf( stderr, "e_init     = %32.16e\n", e_init );

  /* preparation of initial values */
  y_init[1] = p_init;
  y_init[2] = e_init;
  y_init[3] = 0.;
  p[1] = y_init[1];
  e[1] = y_init[2];
  phase[1] = y_init[3];

  /* integration */
  for (count = 2; count <= NUMPOINTS; count++)
    {
      t1 = t_init + (count-2)*Delta;
      t2 = t_init + (count-1)*Delta;
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

  Delta = (t_final - t_init)/(NUMPOINTS-1);
  fprintf(stderr, "delta = %32.16e\n", Delta );

  for (count = 1; count <= NUMPOINTS; count++)
    { 
      t = t_init + (count-1)*Delta;
      fprintf(output_file,"%f %f %f %f %f\n", t,phase[count],p[count],e[count],h[count]);
    }
  fprintf(stderr, "max(t) = %32.16e\n", t);

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







