
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

/*
  convert coefficients of Legendre polynomial on [0, 1] (!!!) to
  coefficients of Bernstein polynomial of the same degree on [0, 1]
  
      Bernstein coef = L2B(degree) %*% Legendre coef

  Farouki (2000) J Comp Appl Math 119 145-160, page 151 (NOT formula 20)
      
*/

SEXP L2B(const SEXP degree) {

    SEXP ans;
    int n, i, j, k, idx;
    double *ret;
    
    n = INTEGER(degree)[0];
    
    PROTECT(ans = allocMatrix(REALSXP, n + 1, n + 1));
    ret = REAL(ans);
    
    for (j = 0; j <= n; j++) {
        for (k = 0; k <= n; k++) {
            idx = k * (n + 1) + j;
            ret[idx] = 0.0;
            for (i = MAX(0, j + k - n); i <= MIN(j, k); i++)
                ret[idx] += pow(-1, k + i) * choose((double) j, (double) i) * 
                                             choose((double) k, (double) i) * 
                                             choose((double) n - j, (double) k - i);
            ret[idx] = ret[idx] / choose((double) n, (double) k);
        }
    }
    UNPROTECT(1);
    return(ans);
}

/*
  convert coefficients of Bernstein polynomial on [0, 1] to
  coefficients of Legendre polynomial of the same degree on [0, 1] (!!!)
  
      Legendre coef = B2L(degree) %*% Bernstein coef

  Farouki (2000) J Comp Appl Math 119 145-160, formula 21
      
*/

SEXP B2L(const SEXP degree) {

    SEXP ans;
    int n, i, j, k, idx;
    double *ret;
    
    n = INTEGER(degree)[0];
    
    PROTECT(ans = allocMatrix(REALSXP, n + 1, n + 1));
    ret = REAL(ans);
    
    for (j = 0; j <= n; j++) {
        for (k = 0; k <= n; k++) {
            idx = k * (n + 1) + j;
            ret[idx] = 0.0;
            for (i = 0; i <= j; i++)
                ret[idx] += pow(-1, i + j) * choose((double) j, (double) i)*choose((double) j, (double) i) / 
                                             choose((double) n + j, (double) k + i);
            ret[idx] = ret[idx] * (2 * j + 1) * choose((double) n, (double) k) / (n + j + 1);
        }
    }
    UNPROTECT(1);
    return(ans);
}

