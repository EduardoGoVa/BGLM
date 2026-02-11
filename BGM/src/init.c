#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP fB_DLN_G(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fB_P(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fBj_DLN_G(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fBj_DLN_G_mtme(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fBj_P(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fBj_P_mtme(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fllp_PLN(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fllp_PLN_multi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fllp_qmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fU_DLN_G(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fUi(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fUi_DLN_G(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fUi_mtme(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP fXb_multi_prod(SEXP, SEXP, SEXP, SEXP);
extern SEXP kron(SEXP, SEXP);
extern SEXP kron_vec(SEXP, SEXP);
extern SEXP prod(SEXP, SEXP);
extern SEXP rtmvnorm_gibbs(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"fB_DLN_G",       (DL_FUNC) &fB_DLN_G,       10},
    {"fB_P",           (DL_FUNC) &fB_P,           10},
    {"fBj_DLN_G",      (DL_FUNC) &fBj_DLN_G,       8},
    {"fBj_DLN_G_mtme", (DL_FUNC) &fBj_DLN_G_mtme,  9},
    {"fBj_P",          (DL_FUNC) &fBj_P,           9},
    {"fBj_P_mtme",     (DL_FUNC) &fBj_P_mtme,     10},
    {"fllp_PLN",       (DL_FUNC) &fllp_PLN,        8},
    {"fllp_PLN_multi", (DL_FUNC) &fllp_PLN_multi,  8},
    {"fllp_qmc",       (DL_FUNC) &fllp_qmc,        7},
    {"fU_DLN_G",       (DL_FUNC) &fU_DLN_G,        9},
    {"fUi",            (DL_FUNC) &fUi,             6},
    {"fUi_DLN_G",      (DL_FUNC) &fUi_DLN_G,       7},
    {"fUi_mtme",       (DL_FUNC) &fUi_mtme,        7},
    {"fXb_multi_prod", (DL_FUNC) &fXb_multi_prod,  4},
    {"kron",           (DL_FUNC) &kron,            2},
    {"kron_vec",       (DL_FUNC) &kron_vec,        2},
    {"prod",           (DL_FUNC) &prod,            2},
    {"rtmvnorm_gibbs", (DL_FUNC) &rtmvnorm_gibbs, 10},
    {NULL, NULL, 0}
};

void R_init_BGM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}