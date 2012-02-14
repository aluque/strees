#include "complex.h"

#define MIN(_A, _B)   ((_A) < (_B) ? (_A) : (_B))

#define MP_MAXL   128        /* Maximum order of the multipole expansion */

#define FOURPI           12.566370614359172954
#define INVFOURPI        0.079577471545947667884
#define SQRTFOURPI       3.5449077018110320546
#define SQRTINVFOURPI    0.28209479177387814347

#define DIMENSIONS         3     /* Sorry, no fourth dimension yet */

#define X 0
#define Y 1
#define Z 2

/*  These defines are intended to provide a general array interface compatible
    source-compatible with the PyArray interface.  */
#ifdef MP_USE_NUMARRAY
# include "Python.h"
# include "math.h"
# include "Numeric/arrayobject.h"
# define mp_array PyArrayObject
# define mp_array_from_dims PyArray_FromDims
# define mp_CDOUBLE PyArray_CDOUBLE
#endif

extern double prefact[MP_MAXL][MP_MAXL];
extern double trans[MP_MAXL][MP_MAXL];

/* multipol.c */
void mp_calc_prefacts ();
void mp_lplm (int l, int m, double x);
void mp_phase_factor(int m, double x, double y);
mp_array *mp_expand (int maxl, double x0, double y0, double z0, 
		     double dx, double dy, double dz, 
		     int inout, mp_array *array);
double mp_eval (int maxl, mp_array *lm, 
		double x, double y, double z, int inout);
mp_array *mp_eval_array (int maxl, mp_array* lm, 
			 double x0, double y0, double z0, 
			 double dx, double dy, double dz, int inout, double r0,
			 mp_array *array);
mp_array *mp_shift (int maxl, double x, double y, double z, 
		    int inout, mp_array *lm);

