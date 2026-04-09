// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"


// [[register]]
SEXP glm_mcmc(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
	      SEXP Rprobinit, SEXP Rmodeldim,
	      SEXP modelprior,  SEXP betaprior,
		  SEXP positions, SEXP levels, SEXP costs, 
		  SEXP Rbestmodel,  SEXP plocal,
	      SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP Rthin, 
	      SEXP family, SEXP Rcontrol, SEXP Rlaplace, SEXP Rparents
			  )
{
  double *probs, MH=0.0, prior_m=1.0, logmargy, postold, postnew; // shrinkage_m
  int i, m, n, pmodel_old, *bestmodel;
  int mcurrent, n_sure;
  int *counts;
  
  glmstptr *glmfamily;
  glmfamily = make_glmfamily_structure(family);
  
  betapriorptr *betapriorfamily;
  betapriorfamily = make_betaprior_structure(betaprior, family);
  
	int nProtected = 0;
	int nModels=LENGTH(Rmodeldim);
	
    /* ----------------------------------------------------------------
    * NEW: Convert R objects -> GSL
    * ---------------------------------------------------------------*/
	gsl_matrix *POS = sexp_to_gsl_matrix (positions);      /* >>> GSL-ADD */
	gsl_vector *LVL = sexp_to_gsl_vector (levels);         /* >>> GSL-ADD */
	gsl_matrix *costs_mat = sexp_to_gsl_matrix (costs);         /* >>> GSL-ADD */
	int nofvars = LENGTH (levels);
	int n_obs = LENGTH (Y); // Number of observations

	// Rprintf("nModels to sample is %d\n", nModels);
	
	/*
	SEXP ANS = PROTECT(allocVector(VECSXP, 17)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 17)); ++nProtected;
	SEXP Rprobs = PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP MCMCprobs= PROTECT(duplicate(Rprobinit)); ++nProtected;
	SEXP R2 = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP shrinkage = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelspace = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP modeldim =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP counts =  PROTECT(duplicate(Rmodeldim)); ++nProtected;
	SEXP beta = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP se = PROTECT(allocVector(VECSXP, nModels)); ++nProtected;
	SEXP deviance = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP modelprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP priorprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP logmarg = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP sampleprobs = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Q = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;
	SEXP Rintercept = PROTECT(allocVector(REALSXP, nModels)); ++nProtected;

	SEXP NumUnique = PROTECT(allocVector(INTSXP, 1)); ++nProtected;
*/
	
	SEXP ANS = PROTECT(allocVector(VECSXP, 17)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 17)); ++nProtected;
	
	SEXP Rprobs = duplicate(Rprobinit);  /* Initializing PIPs (dummy level) */
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));
	
	SEXP modelspace = allocVector(VECSXP, nModels); 
	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));
	
	/* initialize variables */ 
	
	SEXP logmarg = allocVector(REALSXP, nModels); 
	memset(REAL(logmarg), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 2, logmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));
	
	SEXP modelprobs = allocVector(REALSXP, nModels);  
	memset(REAL(modelprobs), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));
	
	SEXP priorprobs = allocVector(REALSXP, nModels); 
	memset(REAL(priorprobs), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));
	
	SEXP sampleprobs = allocVector(REALSXP, nModels); 
	memset(REAL(sampleprobs), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 5, sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));
	
	SEXP deviance = allocVector(REALSXP, nModels); 
	memset(REAL(deviance), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 6, deviance);
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));
	
	SEXP beta = allocVector(VECSXP, nModels); 
	SET_VECTOR_ELT(ANS, 7, beta);
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));
	
	SEXP se = allocVector(VECSXP, nModels);
	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));
	
	SEXP shrinkage = allocVector(REALSXP, nModels); 
  memset(REAL(shrinkage), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));
	
	SEXP modeldim =  allocVector(INTSXP, nModels); 
	memset(INTEGER(modeldim), 0, nModels * sizeof(int));
	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));
	
	SEXP R2 = allocVector(REALSXP, nModels); 
	memset(REAL(R2), 0.0, sizeof(double) *nModels);
	SET_VECTOR_ELT(ANS, 11, R2);
	SET_STRING_ELT(ANS_names, 11, mkChar("R2"));
	
	SEXP Rcounts =  allocVector(INTSXP, nModels); 
	counts = INTEGER(Rcounts);
	memset(counts, 0, nModels * sizeof(int));
	SET_VECTOR_ELT(ANS, 12, Rcounts);
	SET_STRING_ELT(ANS_names, 12, mkChar("freq"));
	
	SEXP MCMCprobs= duplicate(Rprobinit);
	SET_VECTOR_ELT(ANS, 13, MCMCprobs);
	SET_STRING_ELT(ANS_names, 13, mkChar("probne0.MCMC"));
	
	SEXP NumUnique = allocVector(INTSXP, 1); 
	SET_VECTOR_ELT(ANS, 14, NumUnique);
	SET_STRING_ELT(ANS_names, 14, mkChar("n.Unique"));
	
	SEXP Q = allocVector(REALSXP, nModels); 
	SET_VECTOR_ELT(ANS, 15, Q);
	SET_STRING_ELT(ANS_names, 15, mkChar("Q"));
	
	SEXP Rintercept = allocVector(REALSXP, nModels); 
	SET_VECTOR_ELT(ANS, 16, Rintercept);
	SET_STRING_ELT(ANS_names, 16, mkChar("intercept"));
	
	
	
	setAttrib(ANS, R_NamesSymbol, ANS_names);
	
	


	//get dimsensions of all variables
	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	
	int thin = INTEGER(Rthin)[0];

	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs); /* PIPs pointer */
	n = sortvars(vars, probs, p);
	for (i =n; i <p; i++) REAL(MCMCprobs)[vars[i].index] = probs[vars[i].index];
	for (i =0; i <n; i++) REAL(MCMCprobs)[vars[i].index] = 0.0;
	// int noInclusionIs1 = no_prior_inclusion_is_1(p, probs);

	// fill in the sure things
	int *model = ivecalloc(p);
	for (i = n, n_sure = 0; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
		if (model[vars[i].index] == 1) ++n_sure;
	}

	GetRNGstate();

	NODEPTR tree, branch;
	tree = make_node(-1.0);
	//  Rprintf("For m=0, Initialize Tree with initial Model\n");

	m = 0;
	bestmodel = INTEGER(Rbestmodel);
	INTEGER(modeldim)[m] = n_sure;

	// Rprintf("Create Tree\n");
	branch = tree;
	CreateTree(branch, vars, bestmodel, model, n, m, modeldim, Rparents);
	int pmodel = INTEGER(modeldim)[m];
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
	GetModel_m(Rmodel_m, model, p);
	
	// Initial model fit
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
						Rcontrol, Rlaplace, betapriorfamily,
						positions, levels));

	prior_m  = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);

	// Evaluate logmargy
	logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	// shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),"shrinkage"))[0];

	SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q,Rintercept, 
		prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	UNPROTECT(2);

	// MCMC sampling loop
	int nUnique=0, newmodel=0;
	double *real_model = vecalloc(n);
	int *modelold = ivecalloc(p);
	int old_loc = 0;
	int new_loc;
	pmodel_old = pmodel;
	nUnique=1;
	INTEGER(Rcounts)[0] = 1;
	postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
	memcpy(modelold, model, sizeof(int)*p);
	m = 0;
	int *varin= ivecalloc(p);
	int *varout= ivecalloc(p);
	double problocal = REAL(plocal)[0];
	
	while (nUnique < nModels && m < (INTEGER(MCMC_Iterations)[0] + INTEGER(BURNIN_Iterations)[0])) {
		memcpy(model, modelold, sizeof(int)*p);
		pmodel =  n_sure;

		MH = GetNextModelCandidate(pmodel_old, n, n_sure, model, vars, problocal,
                             varin, varout, Rparents);

		branch = tree;
		newmodel= 0;
		for (i = 0; i< n; i++) {
			int bit =  model[vars[i].index];
			if (bit == 1) {
				if (branch->one != NULL) branch = branch->one;
				else newmodel = 1;
			} else {
				if (branch->zero != NULL)  branch = branch->zero;
				else newmodel = 1;
			}
			pmodel  += bit;
		}

		if (pmodel  == n_sure || pmodel == n + n_sure) {
			MH = 1.0/(1.0 - problocal);
		}
		if (newmodel == 1) {
			new_loc = nUnique;
			PROTECT(Rmodel_m = allocVector(INTSXP,pmodel));
			GetModel_m(Rmodel_m, model, p);

			glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
						Rcontrol, Rlaplace, betapriorfamily, positions, levels));

			prior_m  = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);

			// Evaluate logmargy
			logmargy = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
			//shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),"shrinkage"))[0];

			postnew = logmargy + log(prior_m);
		} else {
		  new_loc = branch->where;
		  postnew =  REAL(logmarg)[new_loc] + log(REAL(priorprobs)[new_loc]);
		}

		MH *= exp(postnew - postold); /* Posterior Model Odds */

		//    Rprintf("MH new %lf old %lf\n", postnew, postold);
		if (unif_rand() < MH) {
		 if (newmodel == 1)  {

			/* David adjusted this if statement */
			if ((m % thin) == 0 && m >= INTEGER(BURNIN_Iterations)[0])  {
			  new_loc = nUnique;
			  INTEGER(Rcounts)[new_loc] = 0; /* David added this, as in glm_mcmc_grow.c */
			  insert_model_tree(tree, vars, n, model, nUnique);
			  INTEGER(modeldim)[nUnique] = pmodel;
				//Rprintf("model %d: %d variables\n", m, pmodel);

			  SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q, Rintercept,
                     prior_m, sampleprobs, logmarg, shrinkage, priorprobs, nUnique);
			  UNPROTECT(2);
			  ++nUnique;
			}
			else UNPROTECT(2);
		 }
			old_loc = new_loc;
			postold = postnew;
			pmodel_old = pmodel;
			memcpy(modelold, model, sizeof(int)*p);
		 } else  {
			if (newmodel == 1) UNPROTECT(2);
		}
		
		/* David added this if statement */
		if ((m % thin) == 0 && m >= INTEGER(BURNIN_Iterations)[0])
			INTEGER(Rcounts)[old_loc] += 1; 
		
		/* Shouldn´t we include this inside the previous if statement? */
		for (i = 0; i < n; i++) {
			// store in opposite order so nth variable is first
			real_model[n-1-i] = (double) modelold[vars[i].index];
			REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		}
	m++;
	}

	for (i = 0; i < n; i++) {
		REAL(MCMCprobs)[vars[i].index] /= (double) m; /* Or by nsamples? */
	}

	// Compute marginal probabilities
	mcurrent = nUnique;
	//		Rprintf("NumUnique Models Accepted %d \n", nUnique);
	compute_modelprobs(modelprobs, logmarg, priorprobs,mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);
	
	INTEGER(NumUnique)[0] = nUnique;
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_VECTOR_ELT(ANS, 13, MCMCprobs);
	
 // Rprintf("Decreasing nModels %d to number of unique models accepted %d \n", nModels, nUnique);
	if (nUnique < nModels) {
	  SET_VECTOR_ELT(ANS, 1, resizeVector(modelspace, nUnique));
	  SET_VECTOR_ELT(ANS, 2, resizeVector(logmarg, nUnique));
	  SET_VECTOR_ELT(ANS, 3, resizeVector(modelprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 4, resizeVector(priorprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 5, resizeVector(sampleprobs, nUnique));
	  SET_VECTOR_ELT(ANS, 6, resizeVector(deviance, nUnique));
	  SET_VECTOR_ELT(ANS, 7, resizeVector(beta, nUnique));
	  SET_VECTOR_ELT(ANS, 8, resizeVector(se, nUnique));
	  SET_VECTOR_ELT(ANS, 9, resizeVector(shrinkage, nUnique));
	  SET_VECTOR_ELT(ANS, 10, resizeVector(modeldim, nUnique));
	  SET_VECTOR_ELT(ANS, 11, resizeVector(R2, nUnique));
	  SET_VECTOR_ELT(ANS, 12, resizeVector(Rcounts, nUnique));
	  SET_VECTOR_ELT(ANS, 15, resizeVector(Q, nUnique));
	  SET_VECTOR_ELT(ANS, 16, resizeVector(Rintercept, nUnique));
	}	  
	
	PutRNGstate();
	UNPROTECT(nProtected);
	//Rprintf("Return\n");

	gsl_matrix_free (POS);                                  
    gsl_vector_free (LVL);                                   
    gsl_matrix_free (costs_mat);                                 

	return(ANS);
}

