#include "complex.h"

#define MIN(_A, _B)   ((_A) < (_B) ? (_A) : (_B))

#define MP_MAXL   128        /* Maximum order of the multipole expansion */

#define FOURPI           12.566370614359172954
#define INVFOURPI        0.079577471545947667884
#define SQRTFOURPI       3.5449077018110320546
#define SQRTINVFOURPI    0.28209479177387814347

#ifndef MAX
#define MAX(x, y) (((x) > (y)) ? (x) : (y))
#endif

#ifndef MIN
#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#endif

#define DIMENSIONS         3     /* Sorry, no fourth dimension yet */

#define X 0
#define Y 1
#define Z 2

/*  These defines are intended to provide a general array interface
    source-compatible with the PyArray interface.  */
#ifdef MP_USE_NUMARRAY
# include "Python.h"
# include "math.h"
# include "numpy/arrayobject.h"
# define mp_intp npy_intp
# define mp_array PyArrayObject
# define mp_array_from_dims PyArray_SimpleNew
# define mp_CDOUBLE PyArray_CDOUBLE
# define mp_DOUBLE PyArray_DOUBLE
# define mp_array_getptr1 PyArray_GETPTR1 
# define mp_array_getptr2 PyArray_GETPTR2 
# define mp_array_getptr3 PyArray_GETPTR3 
# define mp_array_getptr4 PyArray_GETPTR4 
# define mp_array_dim PyArray_DIM
#endif

/*  These are macros to easily access to the data of an array */
#define AR3(_AR, _DATA, _I, _J, _K)  (_DATA + _AR->strides[X] * (_I)	\
				      + _AR->strides[Y] * (_J) \
				      + _AR->strides[Z] * (_K))
#define DAR3(_AR, _I, _J, _K)  AR3(_AR, _AR->data, _I, _J, _K)



extern double prefact[MP_MAXL][MP_MAXL];
extern double trans[MP_MAXL][MP_MAXL];

/* multipol.c */
void mp_init (void);
void mp_calc_prefacts (void);
void mp_lplm (int l, int m, double x);
void mp_phase_factor(int m, double x, double y);
mp_array *mp_expand (int maxl, mp_array *r, mp_array *q, int inout);
double mp_eval (int maxl, mp_array *lm, 
		double x, double y, double z, int inout);
mp_array *mp_eval_array (mp_array* lm, mp_array *r,  int inout);
mp_array *mp_shift (double x, double y, double z, 
		    int inout, mp_array *lm);
mp_array *mp_direct (mp_array *r, mp_array *q, mp_array *reval, double a);
mp_array *mp_field_direct (mp_array *r, mp_array *q, mp_array *reval, double a);

/* efield.c */
void ef_efield (double dx, double dy, double dz, mp_array *phi, 
		mp_array *efield);

