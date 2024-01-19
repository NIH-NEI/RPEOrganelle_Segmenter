%module rpeutil

// This tells SWIG to treat char ** as a special case
%typemap(in) char ** {
  /* Check if is a list */
  if (PyList_Check($input)) {
    long long sz = PyList_Size($input);
    long long i = 0;
    $1 = (char **) malloc((sz+1)*sizeof(char *));
    // $1 = new char *[sz+1];
    for (i = 0; i < sz; i++) {
      PyObject *o = PyList_GetItem($input, i);
      if (PyString_Check(o)) {
      	$1[i] = PyString_AsString(o);
      } else
      if (PyUnicode_Check(o)) {
        $1[i] = PyBytes_AS_STRING(PyUnicode_AsEncodedString(o, "utf-8", "Error ~"));
      } else {
        PyErr_SetString(PyExc_TypeError, "list must contain strings");
        SWIG_fail;
      }
    }
    $1[i] = 0;
  } else {
    PyErr_SetString(PyExc_TypeError, "not a list");
    SWIG_fail;
  }
}

// This cleans up the char ** array we malloc'd before the function call
%typemap(freearg) char ** {
  free($1);
  // delete [] $1;
}

%include "std_vector.i"
%include "std_string.i"

namespace std {
   %template(vectori) vector<int>;
   %template(vectord) vector<double>;
};

%{
    #define SWIG_FILE_WITH_INIT
    #include "rpeutil.h"
%}

%include "../numpy.i"

%init %{
    import_array();
%}

%apply (unsigned char* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned char *mask, int hm, int wm)}
%apply (unsigned short* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned short *data, int hd, int wd)}

%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char *mask3d, int zm3d, int hm3d, int wm3d)}
%apply (unsigned short* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned short *data3d, int zd3d, int hd3d, int wd3d)}

// %apply (int** ARGOUTVIEWM_ARRAY1, int *DIM1) {(int** flatint, int* isz)}
// %apply (double** ARGOUTVIEWM_ARRAY1, int *DIM1) {(double** flatdouble, int* dsz)}

%apply (int* IN_ARRAY2, int DIM1, int DIM2) {(int *rois, int nrois, int roisz)}
%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char *ptmask, int npts, int hptm, int wptm)}



%include "rpeutil.h"

