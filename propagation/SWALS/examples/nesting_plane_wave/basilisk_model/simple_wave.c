#include "grid/cartesian1D.h"
#include "saint-venant.h"

//
// This basilisk model simulates a linear, plane wave propagating through constant depth water
// In theory it should travel with no change in depth; however, without many points per wavelength
// it will dissipate. 
// I used this code to confirm that the numerical dissipation in SWALS was similar when using the same grid size
// and theta value
//

int main()
{
  L0 = 300000.;
  X0 = -L0/3.0;
  G = 9.8;
  N = 400*3; //Number of points along x

  // Limiter parameter
  // By default it is 1.3 (can check by running with theta unset, and comparing to theta=1.3)
  // See utils.h for definitions
  theta = 1.3; 

  run();
}

event init (i = 0)
{

  // Domain and wave pars
  double d0 = 100.0;
  double a0 = 0.001;
  // Wavelength = 38 grid points
  double k0 = 1.0 / ((L0/N) * 38.0);

  foreach(){
    // Velocity
    u.x[] = (x < (-3.25/k0) ? 0.0 :
             x < ( 3.25/k0) ? (sqrt(G/(d0)) * a0 * cos(2.*pi*x*k0)) :
             0.0);
    // Stage (elevation = 0 by default)
    h[] = d0 + sqrt(d0/G) * u.x[];
  }

  // On the first timestep, print flow state info
  //printf ("i = %d t = %g\n", i, t);
  foreach()
    fprintf (stdout, "%12.8e %12.8e %12.8e %12.8e %12.8e\n", t, x, h[], u.x[], zb[]);
  fprintf (stdout, "\n");

}

// The final time should be the time it takes to propagate 100km = 100*1000/sqrt(G*depth)
// Not sure how to define this in basilisk, so I hard code it

//event output (t = (3194.38)) {
//  //foreach()
//  //  fprintf (stdout, "%12.8e %12.8e %12.8e %12.8e %12.8e\n", t, x, h[], u.x[], zb[]);
//  //fprintf (stdout, "\n");
//}

event end (t = (3194.38)) {
  //printf ("i = %d t = %g\n", i, t);
  foreach()
    fprintf (stdout, "%12.8e %12.8e %12.8e %12.8e %12.8e\n", t, x, h[], u.x[], zb[]);
  fprintf (stdout, "\n");
  fprintf(stdout, "# Theta %e\n", theta);
}
