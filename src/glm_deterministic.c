// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"
#include <gsl/gsl_vector.h> // GSL-ADD 
#include <gsl/gsl_matrix.h> // GSL-ADD 


/* ---------- helpers: turn SEXPs into GSL objects ---------------------  */

/* Converts an R logical / integer matrix (column-major) to a gsl_matrix
   in row-major order.  (for positions)          */
gsl_matrix * sexp_to_gsl_matrix(SEXP m)
{
    if (TYPEOF(m) != REALSXP && TYPEOF(m) != INTSXP)
        error("'positions' and 'var.costs' must be integer or numeric matrices");

    if (!isMatrix(m))
        error("'m' must be a matrix");

    int *dim = INTEGER(getAttrib(m, R_DimSymbol));
    const int nrow = dim[0], ncol = dim[1];

    gsl_matrix *M = gsl_matrix_alloc(nrow, ncol);

    if (TYPEOF(m) == INTSXP) {
        int *src = INTEGER(m);  // R matrix is column-major
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i)
                gsl_matrix_set(M, i, j, (double)src[i + nrow * j]);
    } else {  // REALSXP
        double *src = REAL(m);
        for (int j = 0; j < ncol; ++j)
            for (int i = 0; i < nrow; ++i)
                gsl_matrix_set(M, i, j, src[i + nrow * j]);
    }

    return M;
}

/* Converts an R integer / numeric vector to a gsl_vector (for levels)  */
gsl_vector * sexp_to_gsl_vector (SEXP v) /* >>> GSL-ADD */
{
    if (TYPEOF(v) != INTSXP && TYPEOF(v) != REALSXP)
        error("`levels' must be a numeric/integer vector");

    const int len = LENGTH(v);
    gsl_vector *g = gsl_vector_alloc(len);

    if (TYPEOF(v) == INTSXP)
        for (int i = 0; i < len; ++i)
            gsl_vector_set(g, i, (double)INTEGER(v)[i]);
    else
        for (int i = 0; i < len; ++i)
            gsl_vector_set(g, i, REAL(v)[i]);

    return g;
}


int topk(Bit **models, double *prob, int k, struct Var *vars, int n, int p);
void insert_children(int subset, double *list, double *subsetsum,
					 int *queue, int *queuesize, int *tablesize,
					 int *parent, int *pattern, int *position,
					 int *type, char *bits, int  n);
void do_insert(int child, double *subsetsum, int *queue);
int get_next(double *subsetsum, int *queue, int *queuesize);
void set_bits(char *bits, int subset, int *pattern, int *position, int n);
void print_subset(int subset, int rank, Bit **models, Bit *model,
				  double *subsetsum, int *pattern, int *position,
				  int n, struct Var *vars, int p);
int withprob(double p);

// [[register]]
SEXP glm_deterministic(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
		       SEXP Rprobinit, SEXP Rmodeldim, SEXP modelprior, SEXP betaprior,
		       SEXP positions, SEXP levels, SEXP costs,
			   SEXP family, SEXP Rcontrol, SEXP Rlaplace) {
	
	int nProtected = 0;
	int nModels=LENGTH(Rmodeldim);

	glmstptr * glmfamily;
	glmfamily = make_glmfamily_structure(family);

	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);
	
	/* ----------------------------------------------------------------
     * NEW: Convert R objects -> GSL
     * ---------------------------------------------------------------*/
	 gsl_matrix *POS = sexp_to_gsl_matrix (positions);      /* >>> GSL-ADD */
	 gsl_vector *LVL = sexp_to_gsl_vector (levels);         /* >>> GSL-ADD */
	 gsl_matrix *costs_mat = sexp_to_gsl_matrix (costs);    /* >>> GSL-ADD */
 
	 int nofvars = LENGTH(levels); // Number of competing variables (excluding the intercept)
	 int n_obs = LENGTH(Y); // Number of observations
 

	//  Rprintf("Allocating Space for %d Models\n", nModels) ;
	SEXP ANS = PROTECT(allocVector(VECSXP, 14)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 14)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP deviance = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Q = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Rintercept = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	
	memset(REAL(modelprobs), 0.0, sizeof(double) *nModels);
	memset(REAL(priorprobs), 0.0, sizeof(double) *nModels);
	memset(REAL(shrinkage), 0.0, sizeof(double) *nModels);
	memset(REAL(logmarg), 0.0, sizeof(double) *nModels);
	memset(REAL(sampleprobs), 0.0, sizeof(double) *nModels);
	memset(REAL(R2), 0.0, sizeof(double) *nModels);
	memset(REAL(Q), 0.0, sizeof(double) *nModels);
	memset(REAL(Rintercept), 0.0, sizeof(double) *nModels);
	memset(REAL(deviance), 0.0, sizeof(double) *nModels);
	memset(INTEGER(modeldim), 0, sizeof(int) *nModels);
	
;
	double *probs,shrinkage_m,logmargy;

	//get dimensions of all variables
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1]; // includes the intercept
	int k = LENGTH(modelprobs); // no of models to evaluate (including some with 0 prior probability)

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs);
	int n = sortvars(vars, probs, p);

	Bit **models = cmatalloc(k,p);
	int *model = (int *) R_alloc(p, sizeof(int));
	memset(model, 0, p*sizeof(int));
	
	// Counts how many variables (excluding the intercept) were forced in the model
	// (and have a prior inclusion probability used for sampling purposes equal to 1)
	// int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);
	k = topk(models, probs, k, vars, n, p);

	/* now fit all top k models */
	for (int m=0; m < k; m++) {
		int pmodel = 0;
		double pigamma = 1.0;
		for (int j = 0; j < p; j++) {
			model[j] = (int) models[m][j];
			pmodel += (int) models[m][j];
			pigamma *= (double)((int) models[m][j])*probs[j] +
				(1.0 - (double)((int) models[m][j]))*(1.0 -  probs[j]);
		}

		INTEGER(modeldim)[m] = pmodel;

		SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
		GetModel_m(Rmodel_m, model, p); // active variables indices
		
		// glm_fit: list with 2 lists: "fit" or "lpy"
		//evaluate logmargy and shrinkage
		SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights,
						    glmfamily, Rcontrol, Rlaplace, betapriorfamily, positions, levels));
		
		// Allocate a GSL vector of size p
		gsl_vector *index = gsl_vector_alloc(p);
		// Copy values from model to index
		for (size_t i = 0; i < p; ++i) {
			gsl_vector_set(index, i, (double)model[i]);
		}

		// double prior_m  = compute_prior_probs(model,pmodel,p, modelprior, noInclusionIs1);
		double prior_m  = compute_prior_probs_enumeration(index, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);
		gsl_vector_free(index);

		SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q, Rintercept,
                 prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);

		REAL(sampleprobs)[m] = pigamma;
		
	}

	compute_modelprobs(modelprobs, logmarg, priorprobs, k); // Posterior Model Probabilities
	compute_margprobs_old(models, modelprobs, probs, k, p); // PIPs for each predictor (dummy and not factor level)

	/*    freechmat(models,k); */
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));

	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));

	SET_VECTOR_ELT(ANS, 2, logmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));

	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));

	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));

	SET_VECTOR_ELT(ANS, 5,sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));

	SET_VECTOR_ELT(ANS, 6, deviance);
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));

	SET_VECTOR_ELT(ANS, 7, beta);
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));

	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));

	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));

	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));

	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));

	SET_VECTOR_ELT(ANS, 12, Q);
	SET_STRING_ELT(ANS_names, 12, mkChar("Q"));

	SET_VECTOR_ELT(ANS, 13, Rintercept);
	SET_STRING_ELT(ANS_names, 13, mkChar("intercept"));


	setAttrib(ANS, R_NamesSymbol, ANS_names);
	UNPROTECT(nProtected);

	 /* ---------------------------------------------------------------
     *  Cleaning of the GSL objects before returning to R
     * -------------------------------------------------------------- */
	 gsl_matrix_free (POS);                                   /* >>> GSL-ADD */
	 gsl_vector_free (LVL);                                   /* >>> GSL-ADD */
	 gsl_matrix_free (costs_mat);                             /* >>> GSL-ADD */
 
	return(ANS);

}
