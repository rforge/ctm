
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP B2L(SEXP);
extern SEXP L2B(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"B2L", (DL_FUNC) &B2L, 1},
    {"L2B", (DL_FUNC) &L2B, 1},
    {NULL, NULL, 0}
};

void R_init_basefun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
