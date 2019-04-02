#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(lwoifortran)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(woifortran)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

void R_init_calcWOI(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}


