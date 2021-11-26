#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void findex(void *, void *, void *, void *, void *, void *, void *, void *);
extern void GMLE(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void icens_prodlim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void iindexSRC(void *, void *, void *, void *, void *, void *, void *);
extern void IntIndexSRC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void life_table(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loo_comprisk(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void loo_surv(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void neighborhoodSRC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void neighborsSRC(void *, void *, void *, void *, void *);
extern void predict_individual_survival(void *, void *, void *, void *, void *, void *, void *, void *);
extern void pred_index(void *, void *, void *, void *, void *, void *, void *);
extern void prodlim_multistates(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void prodlimSRC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void sindexSRC(void *, void *, void *, void *, void *, void *);
extern void summary_prodlim(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"findex",                      (DL_FUNC) &findex,                       8},
    {"GMLE",                        (DL_FUNC) &GMLE,                        11},
    {"icens_prodlim",               (DL_FUNC) &icens_prodlim,               20},
    {"iindexSRC",                   (DL_FUNC) &iindexSRC,                    7},
    {"IntIndexSRC",                 (DL_FUNC) &IntIndexSRC,                 10},
    {"life_table",                  (DL_FUNC) &life_table,                  14},
    {"loo_comprisk",                (DL_FUNC) &loo_comprisk,                14},
    {"loo_surv",                    (DL_FUNC) &loo_surv,                    12},
    {"neighborhoodSRC",             (DL_FUNC) &neighborhoodSRC,              8},
    {"neighborsSRC",                (DL_FUNC) &neighborsSRC,                 5},
    {"predict_individual_survival", (DL_FUNC) &predict_individual_survival,  8},
    {"pred_index",                  (DL_FUNC) &pred_index,                   7},
    {"prodlim_multistates",         (DL_FUNC) &prodlim_multistates,         22},
    {"prodlimSRC",                  (DL_FUNC) &prodlimSRC,                  29},
    {"sindexSRC",                   (DL_FUNC) &sindexSRC,                    6},
    {"summary_prodlim",             (DL_FUNC) &summary_prodlim,             12},
    {NULL, NULL, 0}
};

void R_init_prodlim(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
