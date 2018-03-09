#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CALC_SIB_NUMBERS(void *, void *, void *, void *, void *, void *, void *);
extern void FPEDERR_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void RPEDERR_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"CALC_SIB_NUMBERS", (DL_FUNC) &CALC_SIB_NUMBERS,  7},
    {"FPEDERR_R",        (DL_FUNC) &FPEDERR_R,        19},
    {"RPEDERR_R",        (DL_FUNC) &RPEDERR_R,        19},
    {NULL, NULL, 0}
};

void R_init_pedantics(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
