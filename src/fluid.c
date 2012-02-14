/*  The routines for solving the convection-diffusion-reaction equations
 */

#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "complex.h"

#define MP_USE_NUMARRAY

#include "streamer3d.h"

static void build_neighs (double *p, int *strides, 
			  double *sneigh[DIMENSIONS][3]);
static void build_pij (double *p, int *strides, pij[DIMENSIONS][3]);


#define THIRD  0.3333333333333L
#define SIXTH  0.1666666666666L
#define PSI(_THETA)  MAX(0, MIN(1, MIN(_THETA, THIRD + SIXTH * (_THETA))))

/*  Indexes in the neighbours array */
#define N_MINUS_ONE 0
#define N_PLUS_ONE  1

/*  Calculates the time-derivatives of sigma (dens. of electrons) from
    given SIGMA and EFIELD.  The output is stored in DSDT. */
void
fl_dsigma (double *dr, mp_array *sigma, mp_array *efield, mp_array *dsdt)
{
  int i, j, k;
  double *p, *ep, *dsp;

  assert (phi->nd == 3);
  assert (efield->nd == 4);

  for (i = 0; i < efield->dimensions[X]; i++) {
    for (j = 0; j < efield->dimensions[Y]; j++) {
      for (k = 0; k < efield->dimensions[Z]; k++) {
	p = (double *) DAR3 (sigma, i, j, k);
	ep = (double *) DAR3 (efield, i, j, k);
	dsp = (double *) DAR3 (dsdt, i, j, k);

	*dsp = divergence(dr, p, phi->strides, ep, efield->strides);
      }
    }
  }
}


/*  Builds the array of ADVECTION terms in a given dimension, which
    is known by the corresponding STRIDE.  FIELDS contains pointers to 
    the fields in the two boundaries for this dimension.
*/
static void
advection_terms (double *p, int stride, double *fields[2], double fluxes[2])
{
  /* Defined in [JCP?, Montijn, Hundsdorfer, Ebert], after eq. (3.6) */
  double pij[3];

  /* Neigbours of sigma */
  double *pneigh[3];

  /* 2. Find the pij's */
  for (k = 0; k < 3; k++) {
    double *zero, *minus, *plus;
    zero = (double *) ((char *) p + (k - 1) * stride);
    minus = (double *) ((char *) zero - stride);
    plus = (double *) ((char *) zero + stride);
	
    pneigh[k] = zero;
    pij[k] = (*zero - *minus) / (*plus - *zero);
  }

  for (k = 0; k < 2; k++) {
    if (*fields[k] > 0) {
      fluxes[k] = (-*fields[k] * (*pneigh[k + 1]) 
		   + PSI (1 / pij[k + 1]) * (*pneigh[k] - *pneigh[k + 1]));
    } else {
      fluxes[k] = (-*fields[k] * (*pneigh[k]) 
		   + PSI (pij[k]) * (pneigh[k + 1] - *pneigh[k]));
    }
  }
}

/* Calculates the sum of the advection terms for all dimensions */
static double
divergence (double *dr, double *p, int *strides, double *ep, int *estrides)
{
  int dim;
  double *fields[2], fluxes[2];
  double total;

  fields[N_PLUS_ONE] = ep;

  total = 0;
  for (dim = 0; dim < DIMENSIONS; dim++) {
    fields[N_MINUS_ONE] = (double *) ((char *) ep - estrides[dim]);
    advection_terms (p, strides[dim], fields, fluxes);
    total += (fluxes[N_MINUS_ONE] - fluxes[N_PLUS_ONE]) / dr[dim];
  }

  return total;
}
