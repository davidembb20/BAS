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

	// int nModels0 = INTEGER(RnModels)[0];  // initial guess on number of models to return
	// int nModels = nModels0;
	int nModels = INTEGER(RnModels)[0];  // initial guess on number of models to return
	double expand = REAL(Rexpand)[0]; // increase to grow vectors  
	
	int nProtected = 0;
	
   /* ----------------------------------------------------------------
    * NEW: Convert R objects -> GSL
    * ---------------------------------------------------------------*/
	gsl_matrix *POS = sexp_to_gsl_matrix (positions);      
	gsl_vector *LVL = sexp_to_gsl_vector (levels);         
	gsl_matrix *costs_mat = sexp_to_gsl_matrix (costs);    
	int nofvars = LENGTH (levels);
	int n_obs = LENGTH (Y); 
  
	SEXP ANS = PROTECT(allocVector(VECSXP, 5)); ++nProtected; 
	SEXP ANS_names = PROTECT(allocVector(STRSXP, 5)); ++nProtected; 
	
	SEXP modelspace = allocVector(VECSXP, nModels); 
	SET_VECTOR_ELT(ANS, 0, modelspace);
	SET_STRING_ELT(ANS_names, 0, mkChar("which"));
	
	SEXP logmarg = allocVector(REALSXP, nModels); 
	memset(REAL(logmarg), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 1, logmarg);
	SET_STRING_ELT(ANS_names, 1, mkChar("logmarg"));
		
	SEXP priorprobs = allocVector(REALSXP, nModels); 
	memset(REAL(priorprobs), 0, nModels * sizeof(double));
	SET_VECTOR_ELT(ANS, 2, priorprobs);
	SET_STRING_ELT(ANS_names, 2, mkChar("priorprobs"));
	
	SEXP Rcounts =  allocVector(INTSXP, nModels); 
	memset(INTEGER(Rcounts), 0, nModels * sizeof(int));
	SET_VECTOR_ELT(ANS, 3, Rcounts);  
	SET_STRING_ELT(ANS_names, 3, mkChar("freq"));  
			
	SEXP NumUnique = allocVector(INTSXP, 1); 
	SET_VECTOR_ELT(ANS, 4, NumUnique); 
	SET_STRING_ELT(ANS_names, 4, mkChar("n.Unique")); 
	
	setAttrib(ANS, R_NamesSymbol, ANS_names);

	SEXP modeldim =  allocVector(INTSXP, nModels); // Before: nModels
	memset(INTEGER(modeldim), 0, 1 * sizeof(int)); // Before: nModels

	int p = INTEGER(getAttrib(X,R_DimSymbol))[1];
	int burnin = INTEGER(BURNIN_Iterations)[0];
	int thin = INTEGER(Rthin)[0];
	int mcmc_size = (INTEGER(MCMC_Iterations)[0] / thin) + 1; // Rounding up
	
	double *probs, prior_m=1.0, logmarg_m, postold, postnew;
	int i, m, n, *bestmodel;
	int n_sure;

	glmstptr *glmfamily;
	glmfamily = make_glmfamily_structure(family);
	betapriorptr *betapriorfamily;
	betapriorfamily = make_betaprior_structure(betaprior, family);

	SEXP Rprobs = duplicate(Rprobinit); 
	probs = REAL(Rprobs); /* PIPs pointer */
	struct Var *vars = (struct Var *) R_alloc(p, sizeof(struct Var)); // Info about the model variables.
	n = sortvars(vars, probs, p); /* n = p - 1, if initprobs = "Uniform", right */
	
	/* Max. capacity needed: in every iteration, I need to perform n marginal likelihood evaluations
	   When I visit a previously visit model, I´ll use the marginal likelihood value stored in aux_tree. */
	int aux_capacity = nModels * n; /* Previously: one model per variable flip: (burnin + mcmc_size) * (n) */
	// Stuff for the auxiliary tree
	SEXP aux_priorprobs  = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected; 
	SEXP aux_logmarg     = PROTECT(allocVector(REALSXP, aux_capacity)); ++nProtected;
	SEXP aux_modeldim 	 = PROTECT(allocVector(INTSXP,  1)); ++nProtected; // Before: aux_capacity
	SEXP aux_modelspace  = PROTECT(allocVector(VECSXP,  aux_capacity)); ++nProtected;

	// fill in the sure things
	int *model = ivecalloc(p);
	for (i = n, n_sure = 0; i < p; i++)  {
		model[vars[i].index] = (int) vars[i].prob;
		if (model[vars[i].index] == 1) ++n_sure;
	}

	GetRNGstate();
	m = 0;
	bestmodel = INTEGER(Rbestmodel);

	// Rprintf("Create Tree\n");
	NODEPTR tree, branch;
	tree = make_node(-1.0);
	branch = tree;
	INTEGER(modeldim)[m] = n_sure;
	CreateTree(branch, vars, bestmodel, model, n, m, modeldim, Rparents); // Incrementa INTEGER(modeldim)[m] 

	// Rprintf("Create Auxiliary Tree\n");
	NODEPTR aux_tree, aux_branch;
	aux_tree = make_node(-1.0);
	aux_branch = aux_tree;
	INTEGER(aux_modeldim)[m] = n_sure; 
	CreateTree(aux_branch, vars, bestmodel, model, n, m, aux_modeldim, Rparents);
	
	int pmodel = INTEGER(modeldim)[m];
	SEXP Rmodel_m =	PROTECT(allocVector(INTSXP, pmodel));
	GetModel_m(Rmodel_m, model, p);
	PrintModel_m(Rmodel_m, model, p); // Added by David
	int nUnique = 0, nUniqueVisited = 0; /* I might store less unique models than those I visited (burn-in and thinning)*/
	
	// Initial model fit
	SEXP glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
		Rcontrol, Rlaplace, betapriorfamily, positions, levels));

	// Evaluate logmarg_m and shrinkage_m
	logmarg_m= REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
	prior_m  = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);
	postold = logmarg_m + log(prior_m);
	INTEGER(Rcounts)[0] = 1; // m = 0
	nUniqueVisited++; // is now 1

	/* Set in Auxiliary Tree */
	Set_less_Model_glm(glm_fit, Rmodel_m, prior_m, aux_logmarg, aux_modelspace, aux_priorprobs, m);
	UNPROTECT(2);

	// Burn-in Sampling loop
	int *modelold = ivecalloc(p);
	memcpy(modelold, model, sizeof(int)*p); // Copying (initial model) to modelold
	int *perm = ivecalloc(n);
	
	int newmodel = 0, old_loc = 0, new_loc, thin_count = 0;
	int aux_old_loc = 0, aux_new_loc;
	int component, oldcomponent, newcomponent;
	double log_ratio, max_denom, log_denom;
	
	// Burn-In Period (old: current model and new: proposal / alternative model that serves as an auxiliary for the next step model)
	/* As there´s no thinning, it can´t be enourmous, if thin > 1 */
	for (int iter = 1; iter < (burnin + 1); iter++) {
			
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
			// Rprintf("vars[component].index]: %d \n", vars[component].index);
			oldcomponent = model [vars[component].index]; // Was this chosen variable already in the model or not (1 or 0)
			// Rprintf ("Already in the model? Y/N %d \n", oldcomponent);
			model [vars[component].index] = 1 - model [vars[component].index];  /* Auxiliary model (one-bit flip)
			necessary for the computation of the next model */
	
			/*  ── Checking if the model was visited already / belongs to the tree ────   */
			aux_branch = aux_tree;        /* start at the root of the tree                      */
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
				// PrintModel_m(Rmodel_m, model, p);

				glm_fit = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
											   Rcontrol, Rlaplace, betapriorfamily, positions, levels));

				logmarg_m = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];
				prior_m = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);
				postnew = logmarg_m + log (prior_m);
				if (!isfinite(postnew)) {

					Rprintf("logmarg_m: %.17g \n", logmarg_m);
					Rprintf("log(prior_m): %.17g \n", log (prior_m));

					warning(
						"Invalid postnew: %.17g",
						postnew
					);

				}

				// Resize auxiliary vectors if capacity exceeded
				if (nUniqueVisited >= aux_capacity) {
					Rprintf ("Entrei no expand durante o burnin");
					aux_capacity   = (int)(expand * aux_capacity);
					
					aux_priorprobs = resizeVector(aux_priorprobs, aux_capacity);
					aux_logmarg    = resizeVector(aux_logmarg, aux_capacity);
					//aux_modeldim   = resizeVector(aux_modeldim, aux_capacity);
					aux_modelspace = resizeVector(aux_modelspace, aux_capacity);
				}

				/* Insert in aux_tree, even if this model won´t be the next step model... */
				aux_new_loc = nUniqueVisited;
				// Rprintf ("Inserting in aux_tree in location... %d \n", nUniqueVisited);
				insert_model_tree (aux_tree, vars, n, model, nUniqueVisited);
				//INTEGER(aux_modeldim)[nUniqueVisited] = pmodel;

				Set_less_Model_glm(glm_fit, Rmodel_m, prior_m, aux_logmarg, aux_modelspace, aux_priorprobs, nUniqueVisited);
				++nUniqueVisited;

				UNPROTECT(2);

			} 
			else {
				aux_new_loc = aux_branch->where;
				// Rprintf ("aux_new_loc: %d \n", aux_new_loc);
				postnew = REAL(aux_logmarg)[aux_new_loc] + log(REAL(aux_priorprobs)[aux_new_loc]);
				
				/* Devia tentar dar print a estes modelos que já foram encontrados; 
				só para ver se o Gibbs Sampler está a funcionar bem... (aux_modelspace)[aux_new_loc] */
			}

			/* Check Appendix A of "On Sampling Strategies in BVS Problems with Large Model Spaces" from Gonzalo*/
            /* If oldcomponent = 0, we have just postnew in the numerator */
			/* If oldcomponent = 1, we have just postold in the numerator */
			/* In the numerator we must have "a", i.e., the model with gamma_idx = 1 */
			max_denom = fmax(postnew, postold);
			log_denom = max_denom + log(exp(postnew - max_denom) + exp(postold - max_denom));
			log_ratio =	(oldcomponent * (postold - postnew) + postnew) - log_denom;
			
			/* Teria de definir ratio = exp(log_ratio)*/
			if (!isfinite(exp(log_ratio)) || exp(log_ratio) < 0.0 || exp(log_ratio) > 1.0) {

				Rprintf("postold: %.10f\n", postold);
				Rprintf("postnew: %.10f\n", postnew);
				Rprintf("oldcomponent: %d\n", oldcomponent);
				Rprintf("post numerator: %.10f\n", oldcomponent * (postold - postnew) + postnew);
				Rprintf ("log denom: %.10f\n", log_denom);
				Rprintf("max post: %.10f\n", max_denom);
				
				warning(
					"Invalid ratio: %.17g",
					exp(log_ratio)
				);

				// if (ratio < 0.0) ratio = 0.0;
				// if (ratio > 1.0) ratio = 1.0;

				// if (!isfinite(ratio))
				// 	ratio = 0.5;
			}
			// Rprintf ("Ratio: %.10f\n", exp(log_ratio)); // in between 0 and 1
            newcomponent = bernoulli_draw (exp(log_ratio)); // Drawing from the full conditional

			if (newcomponent != oldcomponent) { // Not staying in the current model
				aux_old_loc = aux_new_loc;
				postold = postnew;
				if (!isfinite(postold)) {

					warning(
						"Bad assignment: %.17g",
						postnew
					);

				}

				memcpy (modelold, model, sizeof(int)*p); // Copying an array of integers from model to modelold
			} 
			/* else => staying in the same place at the next step
			No need to restore modelold, cause that is done at the beginning of the for loop */
		}

	}

	m++; /* m = 1 */
	/* Main Loop */
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
				//Rprintf("Before aux_tree in idx Loop: %d\n", idx);

				// Copying an array of integers from modelold to model
				memcpy (model, modelold, sizeof(int)*p);

				/* Next Variable Flip (perm[ind] is always different from intercept_pos) */
				component    = perm [idx]; // Randomly chosen index (not variable)
				// Rprintf("vars[component].index]: %d \n", vars[component].index);
				oldcomponent = model [vars[component].index]; // Was this chosen variable already in the model or not (1 or 0)
				// Rprintf ("Already in the model? Y/N %d \n", oldcomponent);
				model [vars[component].index] = 1 - model [vars[component].index];

				/*  ── Checking if the model was visited already / belongs to the tree ────   */
				aux_branch = aux_tree;        /* start at the root of the tree                      */
				newmodel   = 0;           /* assume the model already exists                    */
				pmodel     = n_sure;
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
					// PrintModel_m(Rmodel_m, model, p);

					glm_fit      = PROTECT(glm_FitModel(X, Y, Rmodel_m, Roffset, Rweights, glmfamily,
														Rcontrol, Rlaplace, betapriorfamily, positions, levels));

					logmarg_m     = REAL(getListElement(getListElement(glm_fit, "lpy"),"lpY"))[0];					
					prior_m      = compute_prior_probs_MCMC (model, p, modelprior, POS, nofvars, LVL, costs_mat, n_obs);
					postnew = logmarg_m + log (prior_m);
					if (!isfinite(postnew)) {

						Rprintf("logmarg_m: %.17g \n", logmarg_m);
						Rprintf("log(prior_m): %.17g \n", log (prior_m));

						warning(
							"Invalid postnew: %.17g",
							postnew
						);

					}
					
					// Resize auxiliary vectors if capacity exceeded
					if (nUniqueVisited >= aux_capacity) {
						Rprintf("Entrei no loop para expandir aux_capacity"); 
						aux_capacity   = (int)(expand * aux_capacity);

						aux_priorprobs = resizeVector(aux_priorprobs, aux_capacity);
						aux_logmarg    = resizeVector(aux_logmarg, aux_capacity);
						//aux_modeldim   = resizeVector(aux_modeldim, aux_capacity);
						aux_modelspace = resizeVector(aux_modelspace, aux_capacity);
					}

					aux_new_loc = nUniqueVisited;
					// Rprintf ("Inserting in aux_tree in location... %d \n", nUniqueVisited);
					insert_model_tree (aux_tree, vars, n, model, nUniqueVisited);
					//INTEGER(aux_modeldim)[nUniqueVisited] = pmodel;	

					Set_less_Model_glm(glm_fit, Rmodel_m, prior_m, aux_logmarg, aux_modelspace, aux_priorprobs, nUniqueVisited);
					++nUniqueVisited;

					UNPROTECT(2);
				} 
				else {
					aux_new_loc = aux_branch->where;
					// Rprintf ("aux_new_loc: %d \n", aux_new_loc);
					postnew = REAL(aux_logmarg)[aux_new_loc] + log(REAL(aux_priorprobs)[aux_new_loc]);
				}

				/* Check Appendix A of "On Sampling Strategies in BVS Problems with Large Model Spaces" from Gonzalo*/
				/* If oldcomponent = 0, we have just postnew in the numerator */
				/* If oldcomponent = 1, we have just postold in the numerator */
				/* In the numerator we must have "a", i.e., the model with gamma_idx = 1 */
				max_denom = fmax(postnew, postold);
				log_denom = max_denom + log(exp(postnew - max_denom) + exp(postold - max_denom));
				// double postnum;
				// if (oldcomponent == 1)
				// 	postnum = postold;
				// else
				// 	postnum = postnew;
				log_ratio =	(oldcomponent * (postold - postnew) + postnew) - log_denom;
				// Rprintf("post numerator: %.10f\n", oldcomponent * (postold - postnew) + postnew);
				// Rprintf("max post: %.10f\n", max_denom);
				// Rprintf("postold: %.10f\n", postold);
				// Rprintf("postnew: %.10f\n", postnew);
				/* Teria de definir ratio = exp(log_ratio)*/
				if (!isfinite(exp(log_ratio)) || exp(log_ratio) < 0.0 || exp(log_ratio) > 1.0) {

					Rprintf("postold: %.10f\n", postold);
					Rprintf("postnew: %.10f\n", postnew);
					Rprintf("oldcomponent: %d\n", oldcomponent);
					Rprintf("post numerator: %.10f\n", oldcomponent * (postold - postnew) + postnew);
					Rprintf ("log denom: %.10f\n", log_denom);
					Rprintf("max post: %.10f\n", max_denom);
				
					// post numerator: -nan
					// log denom: -1680.2620397829
					// max post: -1680.2620397829
					// postold: -inf
					// postnew: -1680.2620397829

					warning(
						"Invalid ratio: %.17g",
						exp(log_ratio)
					);

					// if (ratio < 0.0) ratio = 0.0;
					// if (ratio > 1.0) ratio = 1.0;

					// if (!isfinite(ratio))
					// 	ratio = 0.5;
				}
				// Rprintf ("Ratio: %.10f\n", exp(log_ratio)); // in between 0 and 1
				newcomponent = bernoulli_draw (exp(log_ratio)); // Drawing from the full conditional

				if (newcomponent != oldcomponent) { 
					aux_old_loc = aux_new_loc;
					postold = postnew;
					if (!isfinite(postold)) {

						warning(
							"Bad assignment: %.17g",
							postnew
						);

					}
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
			
			if (nUnique >= nModels && m < mcmc_size){

				Rprintf("Entrei no loop para expandir nModels"); 
				// expand nModels and grow result vectors
				nModels = (int) (expand*nModels); //add checks to ensure it is not above max int
				
				modelspace = resizeVector(modelspace, nModels);
				SET_VECTOR_ELT(ANS, 0, modelspace);
				
				logmarg = resizeVector(logmarg, nModels);
				SET_VECTOR_ELT(ANS, 1, logmarg);
								
				priorprobs = resizeVector(priorprobs, nModels);
				SET_VECTOR_ELT(ANS, 2, priorprobs);
											
				Rcounts = resizeVector(Rcounts, nModels);
				SET_VECTOR_ELT(ANS, 3, Rcounts); 

				//modeldim = resizeVector(modeldim, nModels);
			}

			new_loc = nUnique;
			/* Any model in tree belongs also to aux_tree*/
			insert_model_tree (tree, vars, n, modelold, nUnique);
			// INTEGER(modeldim)[nUnique] = pmodel;

			// PROTECT (Rmodel_m = allocVector(INTSXP, pmodel)); // pmodel is the number of active variables in the model
			// GetModel_m (Rmodel_m, modelold, p); // Fill Rmodel_m with indices of active variables
			// PrintModel_m(Rmodel_m, modelold, p);
			// UNPROTECT(1);

			Set_less_Model_gibbs (nUnique,	
				REAL(aux_logmarg)[aux_old_loc], REAL(aux_priorprobs)[aux_old_loc], 
				logmarg, priorprobs,
				VECTOR_ELT(aux_modelspace, aux_old_loc), modelspace); 
			++nUnique;

		} 
		else {
			new_loc = branch->where;
		} 
		
		old_loc = new_loc;	
		INTEGER (Rcounts)[old_loc] += 1;
		m++;
	}

	INTEGER(NumUnique)[0] = nUnique;	
	Rprintf("NumUnique Models Accepted %d \n", nUnique);

	//	Rprintf("Decreasing nModels %d to number of unique models accepted %d \n", nModels, nUnique);
	if (nUnique < nModels) {
	  SET_VECTOR_ELT(ANS, 0, resizeVector(modelspace, nUnique));
	  SET_VECTOR_ELT(ANS, 1, resizeVector(logmarg, nUnique)); 
	  SET_VECTOR_ELT(ANS, 2, resizeVector(priorprobs, nUnique)); 
	  SET_VECTOR_ELT(ANS, 3, resizeVector(Rcounts, nUnique)); 
	}	  

	PutRNGstate();
	UNPROTECT(nProtected);

    gsl_matrix_free (POS);                                  
    gsl_vector_free (LVL);                                   
    gsl_matrix_free (costs_mat);

	return(ANS);
}

