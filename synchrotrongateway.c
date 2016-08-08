#include <gsl/gsl_sf_synchrotron.h>
#include <math.h>
/**
   Wrapper for c library GSL to use in Fortran.
   Compile with 
   gcc -c synchrotrongateway.c
   to create a file which can then be linked to the fortran file.
   
 **/

/* Returns normalized synchrotron function
   (normalized so integral over all x is 1)
 */
void synchrotrongateway_(double *x, double *res){
  if ( *x >= 0.0 && *x < 500.0)
    {
      *res = 9.0*sqrt(3.0)/(8.0*3.1415926535)*gsl_sf_synchrotron_1(*x);
    }
  else
    {
      *res = 0.0;
    }
  //  printf("Input:%e Output:%e\n",*x,*res);
}

