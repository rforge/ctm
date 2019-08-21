
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP pnormMRS(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"pnormMRS", (DL_FUNC) &pnormMRS, 1},
    {NULL, NULL, 0}
};

void R_init_basefun(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
