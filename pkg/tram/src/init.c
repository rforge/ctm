
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP R_pnormMRS(SEXP);
extern SEXP R_inner(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"R_pnormMRS", (DL_FUNC) &R_pnormMRS, 1},
    {"R_inner", (DL_FUNC) &R_inner, 1},
    {NULL, NULL, 0}
};

void R_init_basefun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
