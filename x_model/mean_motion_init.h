#include <stdio.h>

struct ic_params
	{
		double f_init; 		/* initial GW frequency in Hz   */ 
		double M;           /* total mass in solar masses   */
		double eta;			/* symmetric mass ratio			*/
		double e;           /* initial time eccentricity    */
		int pn_order;
	};

double n_init( double n, void *params );
