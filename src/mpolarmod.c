/* numutils:  Several extensions to the Numerics package */

#define MP_USE_NUMARRAY

#include "stdlib.h"
#include "Python.h"
#include "math.h"
#include "numpy/arrayobject.h"
#include "streamer3d.h"
#include "complex.h"

#define TRUE 1
#define FALSE 0

static double eval_expansion (int maxl, PyArrayObject *lm, 
			      double x, double y, double z, int inout);
static int iter_next (int nd, int *indices, int *shape);
static char* element_dbl (PyArrayObject *ar, int *indxs);

char *invok_name;

static PyObject *
expand (PyObject *self, PyObject *args)
{
  int maxl, inout;

  PyArrayObject *r, *q, *res;

  if (!PyArg_ParseTuple(args, "iO!O!i", 
			&maxl,
			&PyArray_Type, &r, 
			&PyArray_Type, &q, 
			&inout))
    return NULL;

  if (maxl > MP_MAXL) {
     PyErr_SetString(PyExc_ValueError, 
		     "Not enough memory to compute so many terms.  Change MP_MAXL and recompile");
     return NULL;
  }

  res = mp_expand (maxl, r, q, inout);

  return PyArray_Return (res);
}

/* Evaluates a multipolar expansion at a given point. */
static PyObject *
evaluate (PyObject *self, PyObject *args)
{
  int maxl, inout;
  double x, y, z, res;

  PyArrayObject *array;

  if (!PyArg_ParseTuple(args, "i(ddd)iO!", 
			&maxl, &x, &y, &z, 
			&inout,
			&PyArray_Type, &array))
    return NULL;

  if (array->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "The array of coefficients must have dimension 2.");
     return NULL;
  }

  res = (double) mp_eval (maxl, array, x, y, z, inout);

  return PyFloat_FromDouble (res);
}

static PyObject *
evaluate_array (PyObject *self, PyObject *args)
{
  int inout;

  PyArrayObject *r, *lm, *res;

  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &lm, 
			&PyArray_Type, &r,
			&inout))
    return NULL;

  if (lm->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "The array of coefficients must have dimension 2.");
     return NULL;
  }

  res = mp_eval_array (lm, r, inout);
  
  return PyArray_Return (res);
}

static PyObject *
evaluate_field_array (PyObject *self, PyObject *args)
{
  int inout;

  PyArrayObject *r, *lm, *res;

  if (!PyArg_ParseTuple(args, "O!O!i",
			&PyArray_Type, &lm, 
			&PyArray_Type, &r,
			&inout))
    return NULL;

  if (lm->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "The array of coefficients must have dimension 2.");
     return NULL;
  }

  res = mp_eval_field_array (lm, r, inout);
  
  return PyArray_Return (res);
}


/* Shifts a multipolar expansion to a given point. */
static PyObject *
shift_expansion (PyObject *self, PyObject *args)
{
  int inout;
  double x, y, z;

  PyArrayObject *array, *res;

  if (!PyArg_ParseTuple(args, "(ddd)iO!", 
			&x, &y, &z, 
			&inout,
			&PyArray_Type, &array))
    return NULL;

  if (array->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "The array of coefficients must have dimension 2.");
     return NULL;
  }

  res = mp_shift (x, y, z, inout, array);

  return PyArray_Return (res);
}

/* Shifts a multipolar expansion to a given point. */
static PyObject *
direct (PyObject *self, PyObject *args)
{
  PyArrayObject *r, *q, *reval, *phi;
  double a;

  if (!PyArg_ParseTuple(args, "O!O!O!d", 
			&PyArray_Type, &r,
			&PyArray_Type, &q,
			&PyArray_Type, &reval,
			&a))
    return NULL;

  if (r->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "r must have dimension 2.");
     return NULL;
  }

  if (reval->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "reval must have dimension 2.");
     return NULL;
  }

  if (r->dimensions[1] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "r must have shape (N, 3).");
     return NULL;
  }

  if (reval->dimensions[1] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "reval must have shape (N, 3).");
     return NULL;
  }

  phi = mp_direct (r, q, reval, a);

  return PyArray_Return (phi);
}


/* Shifts a multipolar expansion to a given point. */
static PyObject *
field_direct (PyObject *self, PyObject *args)
{
  PyArrayObject *r, *q, *reval, *efield;
  double a;

  if (!PyArg_ParseTuple(args, "O!O!O!d", 
			&PyArray_Type, &r,
			&PyArray_Type, &q,
			&PyArray_Type, &reval,
			&a))
    return NULL;

  if (r->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "r must have dimension 2.");
     return NULL;
  }

  if (reval->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "reval must have dimension 2.");
     return NULL;
  }

  if (r->dimensions[1] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "r must have shape (N, 3).");
     return NULL;
  }

  if (reval->dimensions[1] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "reval must have shape (N, 3).");
     return NULL;
  }

  efield = mp_field_direct (r, q, reval, a);

  return PyArray_Return (efield);
}


/* Checks if two boxes are near-neighbours by looking at their integer
   coordinates. */
static PyObject *
are_near_neighbours (PyObject *self, PyObject *args)
{
  PyArrayObject *c0, *c1;
  int i, *p0, *p1;

  if (!PyArg_ParseTuple(args, "O!O!", 
			&PyArray_Type, &c0,
			&PyArray_Type, &c1))
    return NULL;

  if (c0->nd != 1) {
     PyErr_SetString(PyExc_ValueError, 
		     "c0 must have dimension 1.");
     return NULL;
  }

  if (c1->nd != 1) {
     PyErr_SetString(PyExc_ValueError, 
		     "c1 must have dimension 1.");
     return NULL;
  }

  if (c0->dimensions[0] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "c0 must have shape (3,).");
     return NULL;
  }

  if (c1->dimensions[0] != 3) {
     PyErr_SetString(PyExc_ValueError, 
		     "c1 must have shape (3,).");
     return NULL;
  }


  for (i = 0; i < 3; i++) {
    p0 = (int *) PyArray_GETPTR1 (c0, i);
    p1 = (int *) PyArray_GETPTR1 (c1, i);
    
    if (*p0 - *p1 > 1 || *p1 - *p0 > 1) Py_RETURN_FALSE;
  }

  Py_RETURN_TRUE;
}


/* Adds two (complex) matrices, overwritting the first one. */
static PyObject *
accum (PyObject *self, PyObject *args)
{
  int i, j;
  PyArrayObject *m0, *m1;
  complex double *p0, *p1;

  if (!PyArg_ParseTuple(args, "O!O!", 
			&PyArray_Type, &m0,
			&PyArray_Type, &m1))
    return NULL;

  if (m0->dimensions[0] < m1->dimensions[0] 
      || m0->dimensions[1] < m1->dimensions[1]) {
     PyErr_SetString(PyExc_ValueError, 
		     "accum: m0 must be able to accomodate m1");
     return NULL;
  }

  if (m0->nd != 2 || m1->nd != 2) {
     PyErr_SetString(PyExc_ValueError, 
		     "accum: the arrays must have dimension 2.");
     return NULL;
  }


  for (i = 0; i < m0->dimensions[0]; i++) {
    /* We make use here of the fact that the matrices come from a mp
       expansion where m <= l. */
    for (j = 0; j <= i; j++) {
      p0 = (complex double *) PyArray_GETPTR2 (m0, i, j);
      p1 = (complex double *) PyArray_GETPTR2 (m1, i, j);

      *p0 += *p1;
    }
  }


  Py_RETURN_NONE;
}



/* Obtains the electric field from the electrostatic potential */
static PyObject *
electric_field (PyObject *self, PyObject *args)
{
  double dx, dy, dz;
  PyArrayObject *phi, *efield;

  if (!PyArg_ParseTuple(args, "(ddd)O!O!",
			&dx, &dy, &dz, 
			&PyArray_Type, &phi,
			&PyArray_Type, &efield))
    return NULL;

  
  ef_efield (dx, dy, dz, phi, efield);
  return PyFloat_FromDouble (0.0);
  /*return PyArray_Return (efield);*/
}

static PyMethodDef mpMethods[] = {
    {"expand",  expand, METH_VARARGS,
     "Computes the multipolar expansion of a charge distribution"},
    {"evaluate",  evaluate, METH_VARARGS,
     "Evauates a multipolar expansion at a given point"},
    {"eval_array",  evaluate_array, METH_VARARGS,
     "Evauates a multipolar expansion for all points of a given array"},
    {"eval_field_array",  evaluate_field_array, METH_VARARGS,
     "Evauates a field from a multipolar expansion for all points of a given array"},
    {"shift",  shift_expansion, METH_VARARGS,
     "Shifts a multipolar expansion to a new point"},
    {"direct",  direct, METH_VARARGS,
     "Directly computes the potential created by a set of charges"},
    {"are_near_neighbours",  are_near_neighbours, METH_VARARGS,
     "Checks if two boxes are near neighbours"},
    {"accum",  accum, METH_VARARGS,
     "Adds two complex matrices and stores the result in the first one"},
    {"field_direct",  field_direct, METH_VARARGS,
     "Directly calculates the electric field from a set of charges"},
    {NULL, NULL, 0, NULL}        /* Sentinel */
};


PyMODINIT_FUNC
initmpolar (void)
{
    PyObject *m;

    m = Py_InitModule ("mpolar", mpMethods);
    invok_name = "mpolar";

    import_array ();
    mp_init ();
}

