#include <stdio.h>
#include <math.h>

int main( int argc, char* argv[] )
{
   FILE* afp = NULL;
   const int fname_len = 256;
   char fname[ fname_len];

   int i, i_max;
   double t, h;
   double x, y, z;
   
   snprintf( fname, fname_len * sizeof(char), "analytic_soln.txt" );
   afp = fopen( fname, "w" ); 

   i = 0;
   i_max = 10000;
   h = 0.002;
   t = 0.0;
 
    for (i = 0; i < i_max; i++)
      {
        t = (double)i * h;
        x = (1./4.) * ( 3. -  exp( -4.* t / 5. ) );         
	y = 1. - exp( -4.* t / 5. ) + exp( -t );
 	z = -(1./2.) * ( 1. - (1./3.) * exp( -4.* t / 5. ) );  	
	
	fprintf( afp, "%.5e %.5e %.5e %.5e \n", t, x, y, z );
      }  	
   fclose(afp);
   return 0;
}
