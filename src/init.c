// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>



/* .C calls */
extern void gexpectations_vect(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hypergeometric1F1(void *, void *, void *, void *, void *, void *);
extern void hypergeometric2F1(void *, void *, void *, void *, void *);
extern void logHyperGauss2F1(void *, void *, void *, void *, void *);
extern void phi1(void *, void *, void *, void *, void *, void *, void *, void*, void*);
extern void tcch(void *, void *, void *, void *, void *, void *, void *, void*);

/* .Call calls */
//extern SEXP glm_fit(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

// David: added 3 extra arguments (type SEXP) in the call (positions, levels, costs) and gibbs functions
extern SEXP glm_deterministic(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_mcmc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_mcmc_grow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_mcmcbas(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_sampleworep(SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_sampleworep_grow(SEXP, SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// added new functions
extern SEXP glm_gibbssampler(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_gibbsBVS(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_gibbssampler_grow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP glm_gibbsBVS_grow(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// ---------------------------------------------------------------------------------

extern SEXP deterministic(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
// extern SEXP mcmc(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mcmc_grow(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP amcmc(SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP amcmc_grow(SEXP, SEXP, SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP mcmcbas(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sampleworep_new(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP sampleworep_grow(SEXP,SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);



static const R_CMethodDef CEntries[] = {
  {"gexpectations_vect", (DL_FUNC) &gexpectations_vect, 11},
  {"hypergeometric1F1",  (DL_FUNC) &hypergeometric1F1,  6},
  {"hypergeometric2F1",  (DL_FUNC) &hypergeometric2F1,  5},
  {"logHyperGauss2F1",   (DL_FUNC) &logHyperGauss2F1,   5},
  {"phi1",               (DL_FUNC) &phi1,               9},
  {"tcch",               (DL_FUNC) &tcch,               8},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"glm_deterministic", (DL_FUNC) &glm_deterministic, 15},
  //{"glm_fit",           (DL_FUNC) &glm_fit,            7},
  {"glm_mcmc",          (DL_FUNC) &glm_mcmc,          20},
  {"glm_mcmc_grow",     (DL_FUNC) &glm_mcmc_grow,     21}, 
  {"glm_mcmcbas",       (DL_FUNC) &glm_mcmcbas,       20},
  {"glm_sampleworep",   (DL_FUNC) &glm_sampleworep,   18},
  {"glm_sampleworep_grow",(DL_FUNC) &glm_sampleworep_grow, 19},
  {"glm_gibbssampler",  (DL_FUNC) &glm_gibbssampler,  20},
  {"glm_gibbsBVS",      (DL_FUNC) &glm_gibbsBVS,      19},
  {"glm_gibbssampler_grow",(DL_FUNC) &glm_gibbssampler_grow, 21},
  {"glm_gibbsBVS_grow", (DL_FUNC) &glm_gibbsBVS_grow, 20},
//  {"mcmc",              (DL_FUNC) &mcmc,            20},
  {"mcmc_grow",         (DL_FUNC) &mcmc_grow,         21},
  {"amcmc",             (DL_FUNC) &amcmc,             21},
  {"amcmc_grow",             (DL_FUNC) &amcmc_grow,   22},
  {"deterministic",     (DL_FUNC) &deterministic,     11},
  {"mcmcbas",           (DL_FUNC) &mcmcbas,           21},
  {"sampleworep_new",   (DL_FUNC) &sampleworep_new,   15},
  {"sampleworep_grow",   (DL_FUNC) &sampleworep_grow, 15},
  {NULL, NULL, 0}
};

void R_init_BAS(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}



