#include "stdlib.h"
#include "stdio.h"
#include "math.h"
#include "complex.h"

#define MP_USE_NUMARRAY

#include "streamer3d.h"
#include "misc.h"

double prefact[MP_MAXL][MP_MAXL];
double trans[MP_MAXL][MP_MAXL];

/*  These arrays are allocated in mp_init and used in many calculations.
    It is not a big waste of memory to have them allocated always to contain
    an expansion of order MP_MAXL.
    TODO:  Why dont allocate them statically?
*/
double complex *phases;
double *lp;

static int iter_next (int nd, int *indices, int *shape);
static char *element_dbl (mp_array *ar, int *indxs);
static void shift_outward (int lmax, double r, mp_array *lm, mp_array *res);
static void shift_inward (int lmax, double r, mp_array *lm, mp_array *res);
static void out_to_in (int maxl, double r, mp_array *lm, mp_array *res);
static double complex iklm (int k, int l, int m);


/* Compute the "prefactors" of the unnormalized spherical harmonics, namely

  prefact[l][m] = sqrt((l-m)!/(l+m)!)
*/
void
mp_calc_prefacts ()
{
  int i, l, m;

  /* First, the diagonal terms */
  prefact[0][0] = 1.0;

  for (i = 1; i < MP_MAXL; i++) {
    prefact[i][i] = prefact[i - 1][i - 1] / sqrt((2.0 * i - 1.0) * (2.0 * i));
  }

  /* Now, the off-diagonal terms. */
  for (m = 0; m < MP_MAXL; m++) {
    for (l = m + 1; l < MP_MAXL; l++) {
      prefact[l][m] = prefact[l - 1][m] * sqrt ((double) (l - m) 
						/ (double) (l + m));
    }
  }
}


/* Compute the coefficients for a translation of a multipole expansion,
   according to (13) of Cheng, Greengard, Rokhlin, J. Comp. Phys., 155, 468,
   (1999).  Namely, we calculate

       A[l, m] = (-1)^l / sqrt((l - m)! (l + m)!)

    However, we sometimes need negative m in A[l, m].  However, l - m
    is always positive.  We will use l and p = l - m as indices and hence
    build

      trans[l][p] = A[l, l - p].
*/
void
mp_calc_trans ()
{
  int l, p;


  for (l = 0; l < MP_MAXL; l++) {
    /* First, the p = 0 terms */
    if (l == 0) {
      trans[0][0] = 1.0;
    } else {
      trans[l][0] = -trans[l - 1][0] / sqrt((2.0 * l - 1) * (2.0 * l));
    }
    
    /* Now the p > 0 terms */
    for (p = 1; p <= 2 * l; p++){
      trans[l][p] = trans[l][p - 1] * sqrt((2.0 * l - p + 1) / p);
    }
  }
}

/*  Compute the associated Legendre functions into lp, which is assumed to
    be an array of at least [l][m] doubles. */
void
mp_lplm (int l, int m, double x)
{
  int i, j, ls;
  double xs, xq;

  for (i = 0; i < l; i++) {
    for (j = 0; j < m; j++) {
      lp[i * m + j] = 0;
    }
  }

  lp[0] = 1.0;

  if (l == 0) return;

  if (fabs (x) == 1.0) {
    for (i = 1; i < l; i++)
      lp[i * m] = lp[(i - 1) * m] * x;
    return;
  }

  ls = fabs(x) > 1.0 ? -1: 1;

  xs = ls * (1.0 - x*x);
  xq = sqrt(xs);
  
  for (i = 1; i < m; i++) {
    lp[i * (m + 1)] = -ls * (2.0 * i - 1.0) * xq * lp[(i - 1) * (m + 1)];
  }

  for (i = 0; i < MIN(m, l - 1); i++) {
    lp[(i + 1) * m + i] = x * lp[i * (m + 1)] * (2.0 * i + 1.0);
  }

  for (i = 0; i < m; i++) {
    for (j = i + 2; j < l; j++) {
      lp[i + j * m] = ((2.0 * j - 1.0) * x * lp[i + (j - 1) * m]
			  - (i + j - 1.0) * lp[i + (j - 2) * m]) / (j - i);
    }
  }
}


/*  Computes the phase factors exp(-i j phi) where phi = atan2(x, y) and
    0 < j < m */
void
mp_phase_factor(int m, double x, double y)
{
  int j;
  double r;
  double complex eiphi;

  r = sqrt (x*x + y*y);
  if (r > 0) {
    eiphi = (x + I * y) / r;
  } else {
    eiphi = 1.0;
  }

  phases[0] = 1.0;
  for (j = 1; j < m; j++) {
    phases[j] = eiphi * phases[j - 1];
  }
}


/*  Calculates the multipole moments of a distribution of k charges given
    in ARRAY.
    maxl               order of the expansion
    r                  array with shape [3, k] with the location of the 
                       charges
    q                  array with shape [k] with the value of the charges
    inout              inout > 0 means outward, inout < 0 means inward
                       expansion.
*/
mp_array *
mp_expand (int maxl, mp_array *r, mp_array *q, int inout)
{
  int i;
  mp_intp dims[2];
  int k, l, m;
  double complex *term;
  mp_array *res;

  term = (double complex *) xcalloc (maxl * maxl, sizeof(double complex));

  /* This is the number of charges; it is taken from q, so if r is
     larger, the remaining charges are simply ignored. */
  k = mp_array_dim(q, 0);

  for (i = 0; i < k; i++) {
    double x, y, z, rho;
    double iq, rpow;
    
    x = * (double*) mp_array_getptr2(r, X, i);
    y = * (double*) mp_array_getptr2(r, Y, i);
    z = * (double*) mp_array_getptr2(r, Z, i);

    iq = * (double*) mp_array_getptr1(q, i);

    rho = sqrt(x * x + y * y + z * z);

    if (rho != 0) {
      mp_lplm (maxl, maxl, z / rho);
      mp_phase_factor (maxl, x, -y);
    } else {
      mp_lplm (maxl, maxl, 0.0);
      mp_phase_factor (maxl, 0.0, 0.0);
    }

    /* inout < 0 means inward expansion; inout > 0 means 
       outward expansion. */
    rpow = (inout > 0)? 1.0: (1.0 / rho);

    for (l = 0; l < maxl; l++) {
      for (m = 0; m <= l; m++) {
	term[l * maxl + m] += lp[l * maxl + m] * phases[m] * iq * rpow;
      }
      rpow = (inout > 0)? (rpow * rho): (rpow / rho);
    } 
  }

  dims[0] = dims[1] = maxl;

  /* res = (mp_array*) (mp_array_from_dims (2, dims, mp_CDOUBLE)); */
  res = (mp_array*) PyArray_ZEROS (2, dims, mp_CDOUBLE, 0);

  for (l = 0; l < maxl; l++) {
    for (m = 0; m <= l; m++) {
      double *p;
      p = (double *) mp_array_getptr2(res, l, m);
      term[l * maxl + m] *= prefact[l][m];
      p[0] = creal(term[l * maxl + m]);
      p[1] = cimag(term[l * maxl + m]);
    }
  }

  free (term);

  return res;
}


/*  Evaluates a multipolar expansion at a given point, with coordinates
    x, y, z relative to the expansion center.
*/
double
mp_eval (int maxl, mp_array *lm, 
	 double x, double y, double z, int inout)
{
  int l, m, s;
  double complex *term;
  double r;
  double sum, rpow;

  r = sqrt(x * x + y * y + z * z);
  
  /* Is this COMPLETELY correct? I think one has to include this in 
     the l,m=0,0 term. */
  if (r == 0) return 0.0;
  
  mp_lplm (maxl, maxl, z / r);
  mp_phase_factor (maxl, x, y);
  
  sum = 0.0;

  rpow = (inout > 0)? (1.0 / r): 1.0;
  for (l = 0; l < maxl; l++) {
    double msum;
    complex double *c;

    c = (complex double *) (lm->data + l * lm->strides[0]);
    msum = lp[l * maxl] * creal (*c) * prefact[l][0];
    for (m = 1; m <= l; m++) {
      c = (complex double *) (lm->data + l * lm->strides[0] + 
			      m * lm->strides[1]);

      /* inout < 0 means inward expansion; inout > 0 means 
	 outward expansion. */
      msum += 2 * lp[l * maxl + m] * creal (phases[m] * (*c)) * prefact[l][m];
    }
    sum += msum * rpow;
    rpow = (inout > 0)? (rpow / r): rpow * r;
  }

  return sum;
}


/* Evaluates a multipolar expansion at each point of a given array and
   adds the result to that array.  The parameters are the same as for
   mp_expand.
*/
mp_array*
mp_eval_array (mp_array* lm, mp_array *r,  int inout)
{
  int i;
  int l, m;
  mp_intp k;
  int maxl;
  double complex *term;
  mp_array *phi;

  /* This is the order of the expansion. */
  maxl = mp_array_dim(lm, 0);

  /* This is the number of points. */
  k = mp_array_dim(r, 1);

  phi = (mp_array*) (mp_array_from_dims (1, (npy_intp*) &k, mp_DOUBLE));

  for (i = 0; i < k; i++) {
    double x, y, z, rho;
    double rpow, sum, *ptr;
    
    x = * (double*) mp_array_getptr2(r, X, i);
    y = * (double*) mp_array_getptr2(r, Y, i);
    z = * (double*) mp_array_getptr2(r, Z, i);

    rho = sqrt(x * x + y * y + z * z);

    if (rho != 0) {
      mp_lplm (maxl, maxl, z / rho);
      mp_phase_factor (maxl, x, y);
    } else {
      mp_lplm (maxl, maxl, 0.0);
      mp_phase_factor (maxl, 0.0, 0.0);
    }
    
    sum = 0.0;

    /* inout < 0 means inward expansion; inout > 0 means 
       outward expansion. */
    rpow = (inout > 0)? (1.0 / rho): 1.0;
    for (l = 0; l < maxl; l++) {
      double msum;
      int sign;
      complex double *c;

      c = (complex double *) (lm->data + l * lm->strides[0]);
      msum = lp[l * maxl] * creal (*c) * prefact[l][0];

      for (m = 1, sign = -1; m <= l; m++, sign = -sign) {
	c = (complex double *) (lm->data + l * lm->strides[0] + 
				m * lm->strides[1]);
	msum += (lp[l * maxl + m] * creal (phases[m] * (*c)) 
		 * prefact[l][m]);

	msum += (lp[l * maxl + m] * creal (conj(phases[m]) * conj(*c)) 
		 * prefact[l][m]);

      }
      sum += msum * rpow;
      rpow = (inout > 0)? (rpow / rho): rpow * rho;
    }

    ptr = (double *) mp_array_getptr1(phi, i);
    *ptr = sum;

  }

  return phi;
}

/* Shifts a multipolar expansion to a given point. 
      inout > 0  -> Shifts an outward expansion
      inout < 0  -> Shifts an inward expansion
      inout = 0  -> Transforms an outward expansion into an inward one.
*/
mp_array *
mp_shift (double x, double y, double z, int inout, mp_array *lm)
{
  int l, m, maxl;
  double r;
  mp_array *res;

  maxl = mp_array_dim(lm, 0);

  /*res = (mp_array*) (mp_array_from_dims (2, lm->dimensions, mp_CDOUBLE));*/
  res = (mp_array*) PyArray_ZEROS (2, lm->dimensions, mp_CDOUBLE, 0);

  r = sqrt(x * x + y * y + z * z);

  if (r == 0) {
    /*  If the two expansions are done around the same point, we just copy
	the input.  Should we better return the ver same array?  As it is
	now, it is safer but slightly slower.  Anyhow, this case is not very
	frequent.  */
    for (l = 0; l < maxl; l++) {
      for (m = 0; m <= l; m++) {
	complex double *c, *to;
	c = (complex double *) (lm->data + l * lm->strides[0] + 
				m * lm->strides[1]);
	to = (complex double *) (res->data + l * res->strides[0] + 
				m * res->strides[1]);
	*to = *c;
      }
    }
    return res;
  }
  
  mp_lplm (maxl, maxl, z / r);
  mp_phase_factor (maxl, x, y);

  
  if (inout > 0) {
    shift_outward (maxl, r, lm, res);
  } else if (inout < 0) {
    shift_inward (maxl, r, lm, res);
  } else {
    out_to_in (maxl, r, lm, res);
  }
  
  return res;
}

/*  Does the hard computations for an OUTWARD expansion. 
    TODO:  Probably, the loops here can be optimized in some way, but I just
    dont see it now.
*/
static void
shift_outward (int maxl, double r, mp_array *lm, mp_array *res)
{
  int j, k; 

  for (j = 0; j < maxl; j++) {
    for (k = 0; k <= j; k++) {
      complex double sum, *to;
      double rpow;
      int l, m;
      rpow = 1.0;
      sum = 0;
      for (l = 0; l <= j; l++) {
	/* TODO:  We can optimize this loop by running only over m >= 0
	   and accounting separately the m > 0 and m < 0 terms. */
	for (m = -l; m <= l; m++) {
	  complex double *c, cval, ph, term;
	  if (abs (k - m) <= (j - l)) {
	    c = (complex double *) mp_array_getptr2(lm, j - l, abs(k - m));
	    cval = ((k - m) > 0) ? (*c): ~(*c);
	    ph = (m >= 0)? ~phases[m]: phases[-m];
	    
	    term = (iklm (k, m, k - m) * rpow
		    * trans[l][l - m] * trans[j - l][j - l - k + m]
		    * lp[l * maxl + abs(m)] * ph * prefact[l][abs(m)] * cval);


	    sum += term;
	  }
	}
	rpow *= r;
      }

      to = (complex double *) mp_array_getptr2(res, j, k); 
      *to = sum / trans[j][j - k];
    }
  }
}

/*  Does the hard computations for an INWARD expansion. */
static void
shift_inward (int maxl, double r, mp_array *lm, mp_array *res)
{
  int j, k; 

  for (j = 0; j < maxl; j++) {
    for (k = 0; k <= j; k++) {
      double complex sum, *to;
      double rpow;
      int l, m;

      /* In th expression we have r^(l - j); this is we start with rpow = 1
	 when l = j. */
      rpow = 1.0;
      sum = 0;
      for (l = j; l < maxl; l++) {
	for (m = -l; m <= l; m++) {
	  complex double *c, cval, ph;
	  if (abs (m - k) <= (l - j)) {
	    c = mp_array_getptr2 (lm, l, abs(m));

	    cval = (m > 0) ? (*c): ~(*c);
	    ph = ((m - k) > 0)? phases[m - k]: ~phases[k - m];
	    
	    sum += iklm (m, m - k, k) * rpow
	      * trans[l - j][l - j - m + k] / trans[l][l - m] 
	      * lp[(l - j) * maxl + abs(m - k)] 
	      * ph * prefact[l - j][abs(m - k)] * cval;
	  }
	}

	rpow *= -r;
      }
      to = (complex double *) mp_array_getptr2 (res, j, k); 
      *to = sum * trans[j][j - k];
    }
  }
}


/*  Does the hard computations for the transformation OUTWARD - INWARD. */
static void
old_out_to_in (int maxl, double r, mp_array *lm, mp_array *res)
{
  int j, k; 
  double rpowj;

  rpowj = 1.0 / r;
  for (j = 0; j < maxl; j++) {
    for (k = 0; k <= j; k++) {
      double complex sum, *to;
      double rpow;
      int l, m;

      rpow = rpowj;
      sum = 0;
      for (l = 0; l < maxl; l++) {
	for (m = -l; m <= l; m++) {
	  complex double *c, cval, ph;
	  if ((l + j) < maxl && abs (m - k) <= (l + j)) {
	    c = mp_array_getptr2 (lm, l, abs(m));

	    cval = (m > 0) ? (*c): ~(*c);
	    ph = ((m - k) > 0)? phases[m - k]: ~phases[k - m];
	    
	    sum += iklm (k - m, k, m) * rpow
	      * trans[l][l - m] / trans[l + j][l + j - m + k] 
	      * lp[(l + j) * maxl + abs(m - k)] 
	      * ph * prefact[l + j][abs(m - k)] * cval;
	  }
	}

	rpow /= -r;
      }
      to = (complex double *) mp_array_getptr2 (res, j, k); 
      *to = sum * trans[j][j - k];
    }
    rpowj /= r;
  }
}


/*  Does the hard computations for the transformation OUTWARD - INWARD. */
static void
out_to_in (int maxl, double r, mp_array *lm, mp_array *res)
{
  int j, k; 
  double rpowj;

  rpowj = 1.0 / r;
  for (j = 0; j < maxl; j++) {
    for (k = 0; k <= j; k++) {
      double complex sum, *to;
      double rpow;
      int l, m;

      rpow = rpowj;
      sum = 0;
      for (l = 0; l < maxl - j; l++) {
	for (m = 0; m <= l; m++) {
	  complex double *c, cval, ph;
	  c = mp_array_getptr2 (lm, l, m);
	  if (abs (m - k) <= (l + j)) {
	    ph = ((m - k) > 0)? phases[m - k]: ~phases[k - m];
	    sum += iklm (k - m, k, m) * rpow
	      * trans[l][l - m] / trans[l + j][l + j - m + k] 
	      * lp[(l + j) * maxl + abs(m - k)] 
	      * ph * prefact[l + j][abs(m - k)] * (*c);
	  }

	  if (m > 0 && (m + k) <= (l + j)) {
	    ph = ~phases[k + m];
	    sum += iklm (k + m, k, m) * rpow
	      * trans[l][l + m] / trans[l + j][l + j + m + k] 
	      * lp[(l + j) * maxl + (m + k)] 
	      * ph * prefact[l + j][m + k] * ~(*c);
	  }
	}

	rpow /= -r;
      }
      to = (complex double *) mp_array_getptr2 (res, j, k); 
      *to = sum * trans[j][j - k];
    }
    rpowj /= r;
  }
}

/* Calculates directly the interaction between two sets of points:
 the sources (r) and the reval.  
 The potential created by a charge q at a distance r is calculated as
 q / (r + a).  This is to allow for "thin-wire" approximations.
 In the multipolar context, note that a is not included in the multipolar
 expansions; hence it is valid only if the smallest multipolar interaction
 is much larger than a.
*/
mp_array *
mp_direct (mp_array *r, mp_array *q, mp_array *reval, double a)
{
  int i, j;
  mp_intp n, m;
  double q0, x, y, z, dx, dy, dz, *phi_ptr;
  mp_array *phi;

  /* This is the number of charges; it is taken from q, so if r is
     larger, the remaining charges are simply ignored. */
  n = mp_array_dim (q, 0);

  /* The number of evaluation points: */
  m = mp_array_dim (reval, 1);

  phi = (mp_array*) PyArray_ZEROS (1, &m, mp_DOUBLE, 0);

  for (j = 0; j < m; j++) {
    x = * (double*) mp_array_getptr2 (reval, X, j);
    y = * (double*) mp_array_getptr2 (reval, Y, j);
    z = * (double*) mp_array_getptr2 (reval, Z, j);
    phi_ptr = (double*) mp_array_getptr1 (phi, j);

    *phi_ptr = 0;

    for (i = 0; i < n; i++) {
      double rn;

      dx = * (double*) mp_array_getptr2 (r, X, i) - x;
      dy = * (double*) mp_array_getptr2 (r, Y, i) - y;
      dz = * (double*) mp_array_getptr2 (r, Z, i) - z;
      q0 = * (double*) mp_array_getptr1 (q, i);
     
      rn = sqrt(dx * dx + dy * dy + dz * dz);
      /* This is important since often the source and measuring points will be
	 the same and we do not like infinities around here. */
      if (a > 0 || rn > 0.0) {
	*phi_ptr += q0 / (rn + a);
      }
    }
  }

  return phi;
}


/* Calculates directly the electric field at reval created by 
   sources q located at r.
*/
mp_array *
mp_field_direct (mp_array *r, mp_array *q, mp_array *reval, double a)
{
  int i, j;
  mp_intp n, m, dim[2];
  double q0, x, y, z, dx, dy, dz, *f_ptr_x, *f_ptr_y, *f_ptr_z;
  mp_array *field;

  /* This is the number of charges; it is taken from q, so if r is
     larger, the remaining charges are simply ignored. */
  n = mp_array_dim (q, 0);

  /* The number of evaluation points: */
  m = mp_array_dim (reval, 1);

  dim[0] = 3;
  dim[1] = m;

  field = (mp_array*) PyArray_ZEROS (2, dim, mp_DOUBLE, 0);

  for (j = 0; j < m; j++) {
    x = * (double*) mp_array_getptr2 (reval, X, j);
    y = * (double*) mp_array_getptr2 (reval, Y, j);
    z = * (double*) mp_array_getptr2 (reval, Z, j);

    f_ptr_x = (double*) mp_array_getptr2 (field, X, j);
    f_ptr_y = (double*) mp_array_getptr2 (field, Y, j);
    f_ptr_z = (double*) mp_array_getptr2 (field, Z, j);

    *f_ptr_x = *f_ptr_y = *f_ptr_z = 0;

    for (i = 0; i < n; i++) {
      double rn;

      dx = * (double*) mp_array_getptr2 (r, X, i) - x;
      dy = * (double*) mp_array_getptr2 (r, Y, i) - y;
      dz = * (double*) mp_array_getptr2 (r, Z, i) - z;
      q0 = * (double*) mp_array_getptr1 (q, i);
     
      rn = sqrt(dx * dx + dy * dy + dz * dz);
      
      /* This is important since often the source and measuring points will be
	 the same and we do not like infinities around here. */
      if (a > 0 || rn > 0.0) {
	*f_ptr_x += dx * q0 / (rn * (rn + a) * (rn + a));
	*f_ptr_y += dy * q0 / (rn * (rn + a) * (rn + a));
	*f_ptr_z += dz * q0 / (rn * (rn + a) * (rn + a));
      }
    }
  }

  return field;
}


/*  This is the function that has to be called from outside to initialize 
    the library.   All initialization code comes here. */
void
mp_init ()
{
  mp_calc_prefacts ();
  mp_calc_trans ();

  lp = (double *) xmalloc (sizeof(double) * MP_MAXL * MP_MAXL);
  phases = (double complex *) xmalloc (sizeof(double complex) * MP_MAXL);

#ifdef MP_USE_NUMARRAY
  import_array ();
#endif /* MP_USE_NUMARRAY */
}


/*  Finaliztion code.  Frees stuff. */
void
mp_finish ()
{
  free (lp);
  free (phases);
}


/*** Helping functions.  They cannot be called from outside this module ***/
/**************************************************************************/

/*  Computes I^(|k| - |l| - |m|) */
static double complex
iklm (int k, int l, int m)
{
  double complex in[4] = {1, I, -1, -I};
  return in[(abs (k) - abs (l) - abs (m)) & 0x3];
}


/*  Iterates a list of INDICES on every element of a matrix with given SHAPE
    returns 1 if there are elements left, 0 if not */
static int
iter_next (int nd, int *indices, int *shape)
{
  int *p, *ps;
  for (p = indices + nd - 1, ps = shape + nd - 1; p >= indices; p--, ps--) {
    if (*p < (*ps) - 1) {
      (*p) ++; 
      break;
    }

    *p = 0;
  }

  if (*indices == 0 && p == indices - 1) {
    return 0;
  } else {
    return 1;
  }
}


/* Returns a pointer to the element of array AR located where specified
   by INDXS */
static char*
element_dbl (mp_array *ar, int *indxs)
{
  char *p;
  int i;

  p = ar->data;

  for (i = 0; i < ar->nd; i++) {
    p += indxs[i] * ar->strides[i];
  }

  return p;
}


#define L 5
#define M 5
int
main ()
{
  double lp[L][M];
  int i, j;

  mp_init ();
  mp_lplm (M, L, 1.0);

  for (i = 0; i < L; i++) {
    for (j = 0; j < M; j++) {
      printf ("%f ", trans[i][j]);
    }
    printf ("\n");
  }

}
