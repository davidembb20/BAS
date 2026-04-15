// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"


// [[register]]
SEXP glm_gibbsBVS_grow(SEXP Y, SEXP X, SEXP Roffset, SEXP Rweights,
	      SEXP Rprobinit, SEXP RnModels, SEXP modelprior, SEXP betaprior,
		  SEXP positions, SEXP levels, SEXP costs, SEXP Rbestmodel, 
		  SEXP BURNIN_Iterations, SEXP MCMC_Iterations, SEXP Rthin, 
	      SEXP family, SEXP Rcontrol, SEXP Rlaplace, SEXP Rparents, SEXP Rexpand)
{

	int nModels0 = INTEGER(RnModels)[0];  // initial guess on number of models to return
	int nModels = nModels0;
	
//	Rprintf("GibbsBVS GROW nModels is %d\n", nModels);
	
	int nProtected = 0;
	int *counts;
	
	double expand = REAL(Rexpand)[0]; // increase to grow vectors  
//  Rprintf("expand is %f\n", expand);

   /* ----------------------------------------------------------------
    * NEW: Convert R objects -> GSL
    * ---------------------------------------------------------------*/
	gsl_matrix *POS = sexp_to_gsl_matrix (positions);      
	gsl_vector *LVL = sexp_to_gsl_vector (levels);         
	gsl_matrix *costs_mat = sexp_to_gsl_matrix (costs);    
	int nofvars = LENGTH (levels);
	int n_obs = LENGTH (Y); 
  
	SEXP ANS = PROTECT(allocVector(VECSXP, 17)); ++nProtected;
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 17)); ++nProtected;
	
	SEXP Rprobs = duplicate(Rprobinit); 
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_STRING_ELT(ANS_names, 0, mkChar("probne0"));
	
	SEXP modelspace = allocVector(VECSXP, nModels); 
	SET_VECTOR_ELT(ANS, 1, modelspace);
	SET_STRING_ELT(ANS_names, 1, mkChar("which"));
	
	SEXP logmarg = allocVector(REALSXP, nModels); 
	memset(REAL(logmarg), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 2, logmarg);
	SET_STRING_ELT(ANS_names, 2, mkChar("logmarg"));
	
	SEXP modelprobs = allocVector(REALSXP, nModels);  
	memset(REAL(modelprobs), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 3, modelprobs);
	SET_STRING_ELT(ANS_names, 3, mkChar("postprobs"));
	
	SEXP priorprobs = allocVector(REALSXP, nModels); 
	memset(REAL(priorprobs), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 4, priorprobs);
	SET_STRING_ELT(ANS_names, 4, mkChar("priorprobs"));
	
	SEXP sampleprobs = allocVector(REALSXP, nModels); 
	memset(REAL(sampleprobs), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 5, sampleprobs);
	SET_STRING_ELT(ANS_names, 5, mkChar("sampleprobs"));
	
	SEXP deviance = allocVector(REALSXP, nModels); 
	memset(REAL(deviance), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 6, deviance);
	SET_STRING_ELT(ANS_names, 6, mkChar("deviance"));
	
	SEXP beta = allocVector(VECSXP, nModels); 
	SET_VECTOR_ELT(ANS, 7, beta);
	SET_STRING_ELT(ANS_names, 7, mkChar("mle"));
	
	SEXP se = allocVector(VECSXP, nModels);
	SET_VECTOR_ELT(ANS, 8, se);
	SET_STRING_ELT(ANS_names, 8, mkChar("mle.se"));
	
	SEXP shrinkage = allocVector(REALSXP, nModels); 
	memset(REAL(shrinkage), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 9, shrinkage);
	SET_STRING_ELT(ANS_names, 9, mkChar("shrinkage"));
	
	SEXP modeldim =  allocVector(INTSXP, nModels); 
	memset(INTEGER(modeldim), 0, nModels * sizeof(int));
	SET_VECTOR_ELT(ANS, 10, modeldim);
	SET_STRING_ELT(ANS_names, 10, mkChar("size"));
	
	SEXP R2 = allocVector(REALSXP, nModels); 
	memset(REAL(R2), 0, nModels * sizeof(double));
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
	memset(REAL(Q), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 15, Q);
	SET_STRING_ELT(ANS_names, 15, mkChar("Q"));
	
	SEXP Rintercept = allocVector(REALSXP, nModels); 
	memset(REAL(Rintercept), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 16, Rintercept);
	SET_STRING_ELT(ANS_names, 16, mkChar("intercept"));
	
	
	setAttrib(ANS, R_NamesSymbol, ANS_names);

	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int burnin = INTEGER(BURNIN_Iterations)[0];
	int thin = INTEGER(Rthin)[0];
	int mcmc_size = (INTEGER(MCMC_Iterations)[0] / thin) + 1; // Rounding up
	
	double *probs, prior_m=1.0, logmarg_m, postold, postnew; // shrinkage_m
	int i, m, n, *bestmodel;
	int mcurrent, n_sure;

	glmstptr *glmfamily;
	glmfamily = make_glmfamily_structure(family);

	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);


	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	probs =  REAL(Rprobs); /* PIPs pointer */
	n = sortvars(vars, probs, p); /* n = p - 1, if initprobs = "Uniform", right */

	Rprintf("n: %d\n", n);
	Rprintf("p: %d\n", p);
	
	/* Max. capacity needed: in every iteration, I need to perform n marginal likelihood evaluations
	   When I visit a previously visit model, I´ll use the marginal likelihood value stored in aux_tree. */
	int aux_capacity = (burnin + mcmc_size) * (n); 
	// Stuff for the auxiliary tree
	SEXP aux_shrinkage   = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_priorprobs  = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected; 
	SEXP aux_logmarg     = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_modeldim 	 = PROTECT(allocVector(INTSXP,  aux_capacity)); ++nProtected; 
	SEXP aux_modelspace  = PROTECT(allocVector(VECSXP,  aux_capacity)); ++nProtected;
	SEXP aux_beta 		 = PROTECT(allocVector(VECSXP,  aux_capacity)); ++nProtected;
	SEXP aux_se 		 = PROTECT(allocVector(VECSXP,  aux_capacity)); ++nProtected;
	SEXP aux_R2 		 = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_deviance    = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_Q           = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_Rintercept  = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;

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

	//  Rprintf("For m=0, Initialize Tree with initial Model\n");
	m = 0;
	bestmodel = INTEGER(Rbestmodel);

	// Rprintf("Create Tree\n");
	NODEPTR tree, branch;
	tree = make_node(-1.0);
	branch = tree;
	INTEGER(modeldim)[m] = n_sure;
	CreateTree(branch, vars, bestmodel, model, n, m, modeldim, Rparents);

	// Rprintf("Create Auxiliary Tree\n");
	NODEPTR aux_tree, aux_branch;
	aux_tree = make_node(-1.0);
	aux_branch = aux_tree;
	INTEGER(aux_modeldim)[m] = n_sure; /* Can´t I just use modeldim here? */
	CreateTree(aux_branch, vars, bestmodel, model, n, m, aux_modeldim, Rparents);
	
	int pmodel = INTEGER(modeldim)[m];
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP,pmodel));
	GetModel_m(Rmodel_m, model, p);
	int nUnique = 0, nUniqueVisited = 0; /* I might store less unique models than those I visited (burn-in and thinning)*/
	
	// Initial model fit
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
		Rcontrol, Rlaplace, betapriorfamily, positions, levels));

	// Evaluate logmarg_m and shrinkage_m
	//logmarg_m= REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	SEXP lpy = getListElement(glm_fit, "lpy");
	SEXP lpY = getListElement(lpy, "lpY");
	if (lpy == R_NilValue || lpY == R_NilValue || TYPEOF(lpY) != REALSXP || LENGTH(lpY) < 1) {
		Rprintf("DEBUG glm_gibbsBVS init: invalid lpy/lpY (pmodel=%d, nUniqueVisited=%d)\n", pmodel, nUniqueVisited);
		Rprintf("DEBUG glm_gibbsBVS init: TYPEOF(lpy)=%d TYPEOF(lpY)=%d LENGTH(lpY)=%d\n",
		        TYPEOF(lpy), TYPEOF(lpY), LENGTH(lpY));
		error("glm_gibbsBVS: invalid glm_fit$lpy$lpY in initial model");
	}
	logmarg_m = REAL(lpY)[0];
	//shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"),"shrinkage"))[0];

	prior_m  = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);
	if (!R_finite(prior_m) || prior_m <= 0.0) {
		error("glm_gibbsBVS: invalid initial prior probability");
	}

	/* Set in Output Tree */ 
	SetModel_glm(glm_fit, Rmodel_m, beta, se, modelspace, deviance, R2, Q, Rintercept, 
		prior_m, sampleprobs, logmarg, shrinkage, priorprobs, m);
	/* Set in Auxiliary Tree */
	SetModel_glm(glm_fit, Rmodel_m, aux_beta, aux_se, aux_modelspace, aux_deviance, aux_R2, aux_Q, 
		aux_Rintercept, prior_m, sampleprobs, aux_logmarg, aux_shrinkage, aux_priorprobs, m);
	
	UNPROTECT(2);

    INTEGER(Rcounts)[0] = 1;
	postold =  REAL(logmarg)[m] + log(REAL(priorprobs)[m]);
	nUnique++; nUniqueVisited++; 

	// Burn-in Sampling loop
	int *modelold = ivecalloc(p);
	memcpy(modelold, model, sizeof(int)*p);
	int *perm = ivecalloc(n);
	double *real_model = vecalloc(n);
	
	int newmodel = 0, old_loc = 0, new_loc, thin_count = 0;
	int aux_old_loc = 0, aux_new_loc;
	
	int component = 1, oldcomponent = 1, newcomponent = 1;
	double ratio = 0.0;
	
	// Burn-In Period (old: current model and new: proposal / alternative model that serves as an auxiliary for the next step model)
	/* As there´s no thinning, it can´t be enourmous, if thin > 1 */
	for (int iter = 1; iter < (burnin + 1); iter++) {
			
		/* It might be more easier for newcomers to use p instead of n? */
		/* Refreshing the permutation vector at each iteration … */		
		for (int j = 0; j < n; ++j) {   
			perm[j] = j;  /* 'j' represents a variable inside vars[] */
		}
			
		randperm (perm, n); /* random permutation: can´t shuffle the intercept */

		for (int idx = 0; idx < n; idx++) { // Loop across all possible variables (except the intercept)...
		
			// Copying an array of integers from modelold to model
			memcpy (model, modelold, sizeof(int)*p);

			/* Next Variable Flip (perm[idx] is always different from intercept_pos) */
			component    = perm [idx]; // Randomly chosen variable
			oldcomponent = model [vars[component].index]; // Before: vars[component].index
			model [vars[component].index] = 1 - model [vars[component].index];  /* Auxiliary model (one-bit flip)
			necessary for the computation of the next model */
	
			/*  ── Checking if the model was visited already / belongs to the tree ────   */
			aux_branch   = aux_tree;        /* start at the root of the tree                      */
			newmodel = 0;           
			pmodel = n_sure;
			for (int i = 0; i < n; i++) { /* Was this "auxiliary" model already visited?              */
				int bit = model[vars[i].index];  /* inclusion flag for current predictor  */

				if (bit == 1) {
					/* Want to follow the 'one' child.  If it is missing, this model has
					* never been stored before — record that fact but keep going so that
					* pmodel is still computed correctly.                                */
					if (aux_branch->one != NULL)
						aux_branch = aux_branch->one;    /* descend one level                     */
					else
						newmodel = 1;            /* missing node → new (unseen) model     */
				} else {
					/* Analogous logic for 'zero' child (variable excluded). */
					if (aux_branch->zero != NULL)
						aux_branch = aux_branch->zero;
					else
						newmodel = 1;
				}

				pmodel += bit;                   /* running count of included variables   */
			}

			if (newmodel == 1) {

				PROTECT (Rmodel_m = allocVector(INTSXP, pmodel)); // pmodel is the number of active variables in the model
				GetModel_m (Rmodel_m, model, p); // Fill Rmodel_m with indices of active variables

				glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
											   Rcontrol, Rlaplace, betapriorfamily, positions, levels));

				lpy = getListElement(glm_fit, "lpy");
				lpY = getListElement(lpy, "lpY");
				if (lpy == R_NilValue || lpY == R_NilValue || TYPEOF(lpY) != REALSXP || LENGTH(lpY) < 1) {
					Rprintf("DEBUG glm_gibbsBVS burnin: invalid lpy/lpY (iter=%d, idx=%d, pmodel=%d, nUniqueVisited=%d)\n",
							iter, idx, pmodel, nUniqueVisited);
					Rprintf("DEBUG glm_gibbsBVS burnin: TYPEOF(lpy)=%d TYPEOF(lpY)=%d LENGTH(lpY)=%d\n",
							TYPEOF(lpy), TYPEOF(lpY), LENGTH(lpY));
					error("glm_gibbsBVS: invalid glm_fit$lpy$lpY during burn-in");
				}
				logmarg_m = REAL(lpY)[0];
				// logmarg_m = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
				// shrinkage_m = REAL(getListElement(getListElement(glm_fit, "lpy"), "shrinkage"))[0];
				prior_m = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);

				postnew = logmarg_m + log (prior_m);
				/* Insert in aux_tree, even if this model won´t be the next step model... */
				aux_new_loc = nUniqueVisited;
				insert_model_tree (aux_tree, vars, n, model, nUniqueVisited);
				INTEGER(aux_modeldim)[nUniqueVisited] = pmodel;
				SetModel_glm(glm_fit, Rmodel_m, aux_beta, aux_se, aux_modelspace, aux_deviance, aux_R2, aux_Q, 
					aux_Rintercept, prior_m, sampleprobs, aux_logmarg, aux_shrinkage, aux_priorprobs, nUniqueVisited);
				++nUniqueVisited;
			
				UNPROTECT(2);

			} 
			else {
				aux_new_loc      = aux_branch->where;
				postnew = REAL(aux_logmarg)[aux_new_loc] + log(REAL(aux_priorprobs)[aux_new_loc]);
			}

			/* Check Appendix A of "On Sampling Strategies in BVS Problems with Large Model Spaces" from Gonzalo*/
            /* If oldcomponent = 0, we have just exp (postnew) in the numerator */
			/* If oldcomponent = 1, we have just exp (postold) in the numerator */
			/* In the numerator we must have "a", i.e., the model with gamma_idx = 1 */
            ratio = (oldcomponent * (exp (postold) - exp (postnew)) + exp (postnew)) / (exp (postnew) + exp (postold));
            newcomponent = bernoulli_draw (ratio); // Drawing from the full conditional

			if (newcomponent != oldcomponent) { // Not staying in the current model

				aux_old_loc = aux_new_loc;
				postold = postnew;
				memcpy (modelold, model, sizeof(int)*p); // Copying an array of integers from model to modelold
			}
			/* else => staying in the same place at the next step
			No need to restore modelold, cause that is done at the beginning of the for loop */
		}

	}


	m++; /* m = 1 */
	while (m < mcmc_size) {
		
		/* Refreshing the permutation vector at each iteration … */		

		for (int j = 0; j < n; ++j) {   
			perm[j] = j;  /* 'j' represents a variable inside vars[] */
		}
		
		/* random permutation of the competing variables: 
			-> doesn´t shuffle the intercept nor other variables forced to be included, as I want all my models with them; 
			this is achieve because we´re working with a different vector... */
		randperm (perm, n); 

		/* Thinning */
		thin_count = 0;
		while (thin_count < thin) {
			
			for (int idx = 0; idx < n; idx++) // Loop across all possible variables...
			{
				// Copying an array of integers from modelold to model
				memcpy (model, modelold, sizeof(int)*p);

				/* Next Variable Flip (perm[ind] is always different from intercept_pos) */
				component    = perm [idx]; // Randomly chosen index (not variable)
				oldcomponent = model [vars[component].index]; // Before: vars[component].index
				model [vars[component].index] = 1 - model [vars[component].index]; 

				/*  ── Checking if the model was visited already / belongs to the tree ────   */
				aux_branch   = aux_tree;        /* start at the root of the tree                      */
				newmodel = 0;           /* assume the model already exists                    */
				pmodel   = n_sure;
				for (int i = 0; i < n; i++) { 
					int bit = model[vars[i].index];  /* inclusion flag for current predictor  */

					if (bit == 1) {
						/* Want to follow the 'one' child.  If it is missing, this model has
						* never been stored before — record that fact but keep going so that
						* pmodel is still computed correctly.                                */
						if (aux_branch->one != NULL)
							aux_branch = aux_branch->one;    /* descend one level                     */
						else
							newmodel = 1;            /* missing node → new (unseen) model     */
					} else {
						/* Analogous logic for 'zero' child (variable excluded). */
						if (aux_branch->zero != NULL)
							aux_branch = aux_branch->zero;
						else
							newmodel = 1;
					}

					pmodel += bit;                   /* running count of included variables   */
				}

				if (newmodel == 1) {

					PROTECT (Rmodel_m = allocVector(INTSXP, pmodel)); // pmodel is the number of active variables in the model
					GetModel_m (Rmodel_m, model, p); // Fill Rmodel_m with indices of active variables

					glm_fit      = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
														Rcontrol, Rlaplace, betapriorfamily, positions, levels));
					lpy = getListElement(glm_fit, "lpy");
					lpY = getListElement(lpy, "lpY");
					if (lpy == R_NilValue || lpY == R_NilValue || TYPEOF(lpY) != REALSXP || LENGTH(lpY) < 1) {
						Rprintf("DEBUG glm_gibbsBVS main: invalid lpy/lpY (m=%d, idx=%d, thin_count=%d, pmodel=%d, nUniqueVisited=%d)\n",
								m, idx, thin_count, pmodel, nUniqueVisited);
						Rprintf("DEBUG glm_gibbsBVS main: TYPEOF(lpy)=%d TYPEOF(lpY)=%d LENGTH(lpY)=%d\n",
								TYPEOF(lpy), TYPEOF(lpY), LENGTH(lpY));
						error("glm_gibbsBVS: invalid glm_fit$lpy$lpY during main loop");
					}
					logmarg_m = REAL(lpY)[0];
					// logmarg_m     = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
					//shrinkage_m  = REAL(getListElement(getListElement(glm_fit, "lpy"), "shrinkage"))[0];
					prior_m      = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);

					postnew = logmarg_m + log (prior_m);

					aux_new_loc = nUniqueVisited;
					insert_model_tree (aux_tree, vars, n, model, nUniqueVisited);
					INTEGER(aux_modeldim)[nUniqueVisited] = pmodel;
					SetModel_glm(glm_fit, Rmodel_m, aux_beta, aux_se, aux_modelspace, aux_deviance, aux_R2, aux_Q, 
						aux_Rintercept, prior_m, sampleprobs, aux_logmarg, aux_shrinkage, aux_priorprobs, nUniqueVisited);
					++nUniqueVisited;
					UNPROTECT(2);
				} 
				else {
					aux_new_loc = aux_branch->where;
					postnew =  REAL(aux_logmarg)[aux_new_loc] + log(REAL(aux_priorprobs)[aux_new_loc]);
				}

				/* Check Appendix A of "On Sampling Strategies in BVS Problems with Large Model Spaces" from Gonzalo*/
				/* If oldcomponent = 0, we have just exp (postnew) in the numerator */
				/* If oldcomponent = 1, we have just exp (postold) in the numerator */
				/* In the numerator we must have "a", i.e., the model with gamma_idx = 1 */
				ratio = (oldcomponent * (exp (postold) - exp (postnew)) + exp (postnew)) / (exp (postnew) + exp (postold));
				newcomponent = bernoulli_draw (ratio); // Drawing from the full conditional

				if (newcomponent != oldcomponent) { 
					aux_old_loc = aux_new_loc;
					postold = postnew;
					memcpy (modelold, model, sizeof(int)*p); // Copying an array of integers from model to modelold
				}

			}

			thin_count++;
		}

		branch   = tree;        /* start at the root of the tree                      */
		newmodel = 0;  
		pmodel   = n_sure;
		for (int i = 0; i < n; i++) { 
			int bit = modelold[vars[i].index];  /* inclusion flag for current predictor  */

			if (bit == 1) {
				/* Want to follow the 'one' child.  If it is missing, this model has
				* never been stored before — record that fact but keep going so that
				* pmodel is still computed correctly.                                */
				if (branch->one != NULL)
					branch = branch->one;    /* descend one level                     */
				else
					newmodel = 1;            /* missing node → new (unseen) model     */
			} else {
				/* Analogous logic for 'zero' child (variable excluded). */
				if (branch->zero != NULL)
					branch = branch->zero;
				else
					newmodel = 1;
			}

			pmodel += bit;                   /* running count of included variables   */
		}

		if (newmodel == 1) {
				
			new_loc = nUnique;
			/* Any model in tree belongs also to aux_tree*/
			insert_model_tree (tree, vars, n, modelold, nUnique);
			INTEGER(modeldim)[nUnique] = pmodel;

			SetModel_gibbs(nUnique,	
				REAL(aux_logmarg)[aux_old_loc], REAL(aux_shrinkage)[aux_old_loc], REAL(aux_priorprobs)[aux_old_loc],
				logmarg, shrinkage, priorprobs, sampleprobs,
				REAL(aux_deviance)[aux_old_loc], REAL(aux_R2)[aux_old_loc], REAL(aux_Q)[aux_old_loc], REAL(aux_Rintercept)[aux_old_loc],
				deviance, R2, Q, Rintercept,
				VECTOR_ELT(aux_beta, aux_old_loc), VECTOR_ELT(aux_se, aux_old_loc), VECTOR_ELT(aux_modelspace, aux_old_loc), 
				beta, se, modelspace); 
			++nUnique;

		} 
		else {
			new_loc = branch->where;
		} 
		
		old_loc = new_loc;
		INTEGER (Rcounts)[old_loc] += 1;

		for (i = 0; i < n; i++) {
			// store in opposite order so nth variable is first
			real_model[n-1-i] = (double) modelold[vars[i].index];
			REAL(MCMCprobs)[vars[i].index] += (double) modelold[vars[i].index];
		}
		
		if (nUnique >= nModels && m < (INTEGER(MCMC_Iterations)[0]/thin)){
		  // expand nModels and grow result vectors
		  nModels = (int) (expand*nModels); //add checks to ensure it is not above max int
		  
//		  Rprintf("Grow vectors:  Number of unique models %d; nModels is now %d\n", nUnique, nModels); // Need to use growable vector here
		  
		  modelspace = resizeVector(modelspace, nModels);
		  SET_VECTOR_ELT(ANS, 1, modelspace);
		  
		  logmarg = resizeVector(logmarg, nModels);
		  SET_VECTOR_ELT(ANS, 2, logmarg);
		  
		  modelprobs = resizeVector(modelprobs, nModels);
		  SET_VECTOR_ELT(ANS, 3, modelprobs);
		  
		  priorprobs = 	resizeVector(priorprobs, nModels);
		  SET_VECTOR_ELT(ANS, 4, priorprobs);
		  
		  sampleprobs = resizeVector(sampleprobs, nModels);
		  SET_VECTOR_ELT(ANS, 5, sampleprobs);
		  
		  deviance = resizeVector(deviance, nModels);
		  SET_VECTOR_ELT(ANS, 6, deviance);
		  
		  beta = resizeVector(beta, nModels);
		  SET_VECTOR_ELT(ANS, 7, beta);
		  
		  se = resizeVector(se, nModels);
		  SET_VECTOR_ELT(ANS, 8, se);
		  
		  shrinkage = resizeVector(shrinkage, nModels);
		  SET_VECTOR_ELT(ANS, 9, shrinkage);
		  
		  modeldim = resizeVector(modeldim, nModels);
		  SET_VECTOR_ELT(ANS, 10, modeldim);
		  
		  R2 = resizeVector(R2, nModels);
		  SET_VECTOR_ELT(ANS, 11, R2);
		  
		  Rcounts = resizeVector(Rcounts, nModels);
		  SET_VECTOR_ELT(ANS, 12, Rcounts);
		  
		  Q = resizeVector(Q, nModels);
		  SET_VECTOR_ELT(ANS, 15, Q);
		  
		  Rintercept = resizeVector(Rintercept, nModels);
		  SET_VECTOR_ELT(ANS, 16, Rintercept);
		}
		
	m++;
	}


	//	Rprintf("Compute MCMC Probabilities\n");
	for (i = 0; i < n; i++) {
		REAL(MCMCprobs)[vars[i].index] /= (double) m;
	}

	// Compute marginal probabilities
	mcurrent = nUnique;
	//		Rprintf("NumUnique Models Accepted %d \n", nUnique);
	compute_modelprobs(modelprobs, logmarg, priorprobs, mcurrent);
	compute_margprobs(modelspace, modeldim, modelprobs, probs, mcurrent, p);

	INTEGER(NumUnique)[0] = nUnique;
	SET_VECTOR_ELT(ANS, 0, Rprobs);
	SET_VECTOR_ELT(ANS, 13, MCMCprobs);
	
	//	Rprintf("Decreasing nModels %d to number of unique models accepted %d \n", nModels, nUnique);
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

    gsl_matrix_free (POS);                                  
    gsl_vector_free (LVL);                                   
    gsl_matrix_free (costs_mat);

	return(ANS);
}

