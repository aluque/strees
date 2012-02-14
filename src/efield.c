/* Routines to calculate the electric field from the electrostatic potential
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "complex.h"

#define MP_USE_NUMARRAY

#include "streamer3d.h"

/* Calculates the components of the electric field from a matrix
   containing the electrostatic potential 

      i-1    i    i+1
    +-----+-----+-----+     The electrostatic potential is defined
    |     |     |     |   in the cell centers (O); the y-component of
j+1 |  O  |  O  |  O  |   the electric field is defined in the cell
    |     |     |     |   boundaries of the y direction (Y) the 
    +-----+--Y--+-----+   x-component is defined in (X) and the z-component
    |     |     |     |   (not shown) is defined in the surface boundaries 
 j  |  O  X  O  X  O  |   perpendicular to the z axis.  Thus, for example
    |     |     |     |
    +-----+--Y--+-----+    E_X[i, j] = (phi[i+1, j] - phi[i, j]) / h
    |     |     |     |
j-1 |  O  |  O  |  O  |
    |     |     |     |
    +-----+-----+-----+

  IMPORTANT: We are assuming that if the shape of phi is phi[I, J, K], the
  shape of efield is efield[I, J, K, 3]; the last index stands for the three
  components of the e-field.
  TODO:  Now the implementation relies in the C ordering of the matrix
  (the last index is the fastest).  One should perhaps not rely on this.

*/
void
ef_efield (double dx, double dy, double dz, mp_array *phi, mp_array *efield)
{
  int i, j, k;
  /* p = phi[i, j, k]; pi = phi[i + 1, j, k]; pj[i, j + 1, k] ... */
  char *p, *ep;

  assert (phi->nd == 3);
  assert (efield->nd == 4);

  p = phi->data;
  ep = efield->data;
    
  for (i = 0; i < efield->dimensions[X]; i++) {
    for (j = 0; j < efield->dimensions[Y]; j++) {
      for (k = 0; k < efield->dimensions[Z]; k++) {
	char *pi, *pj, *pk; 
	if (i != efield->dimensions[X] && j != efield->dimensions[Y]
	    && k != efield->dimensions[Z]) {
	  pi = p + phi->strides[X];
	  pj = p + phi->strides[Y];
	  pk = p + phi->strides[Z];
	
	  ((double*) ep)[X] = -(*(double *) pi - *(double *) p) / dx;
	  ((double*) ep)[Y] = -(*(double *) pj - *(double *) p) / dy;
	  ((double*) ep)[Z] = -(*(double *) pk - *(double *) p) / dz;
	}
      
	ep += efield->strides[Z];
	p += phi->strides[Z];
      }
    }
  }
}
