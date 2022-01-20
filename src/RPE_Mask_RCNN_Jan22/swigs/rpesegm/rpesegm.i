%module rpesegm

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

%{
    #define SWIG_FILE_WITH_INIT
    #include "rpesegm.h"
%}

%include "../numpy.i"

%init %{
    import_array();
%}

%apply (unsigned char* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned char *mask, int hm, int wm)}
%apply (unsigned char* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned char *angle, int ha, int wa)}
%apply (unsigned short* IN_ARRAY2, int DIM1, int DIM2) {(unsigned short *data, int hd, int wd)}
%apply (unsigned char* IN_ARRAY2, int DIM1, int DIM2) {(unsigned char *mask2, int hm2, int wm2)}

%apply (unsigned short* INPLACE_ARRAY2, int DIM1, int DIM2) {(unsigned short *sdata, int shd, int swd)}
%apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {(double *ddata, int dhd, int dwd)}

%apply (unsigned char* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned char *mask3d, int zm3d, int hm3d, int wm3d)}
%apply (unsigned short* INPLACE_ARRAY3, int DIM1, int DIM2, int DIM3) {(unsigned short *data3d, int zd3d, int hd3d, int wd3d)}

%include "rpesegm.h"

%rename (segment_actin) ren_segment_actin;
%rename (segment_dna) ren_segment_dna;
%rename (compare_with_reshape) ren_compare_with_reshape;
%rename (colorize_reshape_comparison) ren_colorize_reshape_comparison;
%rename (segment_actin_z01) ren_segment_actin_z01;

%inline %{
    void ren_segment_actin(
			unsigned char *mask, int hm, int wm,
			unsigned short *data, int hd, int wd,
			unsigned char *angle, int ha, int wa,
			int otsu1, int otsu2) {
    	segment_actin(wd, hd, data, angle, mask, otsu1, otsu2);
	}
    void ren_segment_dna(
			unsigned char *mask, int hm, int wm,
			unsigned short *data, int hd, int wd,
			unsigned char *angle, int ha, int wa,
			int otsu1, int otsu2) {
    	segment_dna(wd, hd, data, angle, mask, otsu1, otsu2);
	}
	void ren_compare_with_reshape(
			unsigned char *mask, int hm, int wm,
			unsigned char *angle, int ha, int wa,
			const char *rs_csvfile,
			const char *out_csvfile) {
		compare_with_reshape(wm, hm, mask, angle, rs_csvfile, out_csvfile);
	}
	void ren_colorize_reshape_comparison(
			unsigned char *mask3d, int zm3d, int hm3d, int wm3d,
			unsigned short *data, int hd, int wd,
			unsigned char *mask, int hm, int wm) {
		colorize_reshape_comparison(wd, hd, mask3d, data, mask);
	}
    void ren_segment_actin_z01(
			unsigned char *mask, int hm, int wm,
			unsigned short *data, int hd, int wd,
			unsigned char *angle, int ha, int wa,
			unsigned char *mask2, int hm2, int wm2,
			int otsu1, int otsu2) {
    	segment_actin_z01(wd, hd, data, angle, mask, mask2, otsu1, otsu2);
	}
%}
