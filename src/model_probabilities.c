// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"

/* Posterior Model Probabilities computation */
void compute_modelprobs(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, int k)
{
	int m;
	double nc, bestmarg, *modelprobs, *logmarg, *priorprobs;

	logmarg = REAL(Rlogmarg);
	modelprobs = REAL(Rmodelprobs); /* Output Vector */
	priorprobs = REAL(Rpriorprobs);

	bestmarg = logmarg[0]; /* Highest log marginal likelihood sampled (not enumeration) or possible (enumeration) */
	nc = 0.0; /* Normalizing Constant */

	for (m = 0; m < k; m++) {
		if (logmarg[m] > bestmarg) bestmarg = logmarg[m]; /* Appears both in the numerator and denominator*/
	}

	for (m = 0; m < k; m++) {
		modelprobs[m] = logmarg[m] - bestmarg; /* Renormalized log marginal likelihoods*/
		nc += exp(modelprobs[m])*priorprobs[m]; /* Denominator*/
	}

  /* "Bayes Theorem"; things are computed in log scale 1st for stability */
	for (m = 0; m < k; m++) {
	  /*		modelprobs[m] = exp(modelprobs[m] +
			log(priorprobs[m]) - log(nc)); */
		modelprobs[m] = exp(modelprobs[m] - log(nc))*priorprobs[m];
	}
}

/* Used in lm_amcmc*.c */
void compute_modelprobs_HT(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, 
                           SEXP Rsampleprobs, int k, int MCsamples)
{
  int m;
  double nc, bestmarg, *modelprobs, *logmarg, *priorprobs, *sampleprobs;
  
  logmarg = REAL(Rlogmarg);
  modelprobs = REAL(Rmodelprobs);
  priorprobs = REAL(Rpriorprobs);
  sampleprobs = REAL(Rsampleprobs); 
  bestmarg = logmarg[0];
  nc = 0.0;
  
  for (m = 0; m < k; m++) {
    if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
    if (sampleprobs[m] > 0.0) modelprobs[m]  =  -log(1.0 - (pow(1.0 - sampleprobs[m], (double) MCsamples)));
  }
  
  for (m = 0; m < k; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] += logmarg[m] - bestmarg;
      nc += exp(modelprobs[m])*priorprobs[m];
    }
  }
  
  for (m = 0; m < k; m++) {
    if (sampleprobs[m] > 0.0) modelprobs[m] = exp(modelprobs[m] - log(nc))*priorprobs[m];
    else {modelprobs[m] = 0.0;}
  }
}

/* Used in lm_mcmcbas.c*/
void compute_modelprobs_Bayes_HT(SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, 
                           SEXP Rsampleprobs, int M, double *eta, double *nc)
{
  int m;
  double bestmarg, *modelprobs, *logmarg, *priorprobs, *sampleprobs, 
         HT = 0.0, probinS = 0.0, correction = 0.0, nmodels = 0.0;
  
  logmarg = REAL(Rlogmarg);
  modelprobs = REAL(Rmodelprobs);
  priorprobs = REAL(Rpriorprobs);
  sampleprobs = REAL(Rsampleprobs); 
  bestmarg = logmarg[0];

  
  for (m = 0; m < M; m++) {
    if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
  }
  
  for (m = 0; m < M; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] = logmarg[m] - bestmarg + log(priorprobs[m]);
      probinS += sampleprobs[m];
      HT += exp(modelprobs[m] - log(sampleprobs[m]));  
      *nc += exp(modelprobs[m]);
      nmodels += 1.0;
    }
  }  
    
  *eta = HT/nmodels;
  correction = *eta*(1.0 - probinS);
  Rprintf("eta = %lf probinS = %lf  NC = %lf  correction = %lf", *eta, probinS, *nc, correction);
  *nc += (1.0 - probinS)* *eta;
  Rprintf(" corrected NC = %lf \n", *nc);
  
  for (m = 0; m < M; m++) {
    if (sampleprobs[m] > 0.0) {
      modelprobs[m] = exp(modelprobs[m] - log(*nc));
    }
    else {modelprobs[m] = 0.0;}
    
  }
}

void compute_margprobs(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, double *margprobs, 
                       int k, int p) {
	int m, j, *model;
	double *modelprobs;
	modelprobs = REAL(Rmodelprobs);
	for (j=0; j< p; j++)  margprobs[j] = 0.0;
	for(m=0; m< k; m++) {
		model = INTEGER(VECTOR_ELT(modelspace,m));
		for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
			margprobs[model[j]] += modelprobs[m];
		}
	}
}

/* Used in lm_mcmcbas.c*/
void compute_margprobs_Bayes_BAS_MCMC(SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, SEXP Rprobs, SEXP Rsampleprobs, 
                       int M, int p, double eta, double NC)
{
  int m, j, *model, *n;
  double *modelprobs;
  double *beta, probNotInS = 1.0;
  SEXP samplemargs = PROTECT(duplicate(Rprobs)); 
  modelprobs = REAL(Rmodelprobs);
  
  for (j=0; j< p; j++) {
   // Rprintf("j = %d uncorrected pip = %lf \n", j, REAL(Rprobs)[j]);
    REAL(Rprobs)[j] = 0.0;

  }
  modelprobs = REAL(Rmodelprobs);
  for(m=0; m < M; m++) {
    model = INTEGER(VECTOR_ELT(modelspace,m));
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
      REAL(Rprobs)[model[j]] += modelprobs[m];
    }
  }
  
  for (j = 0; j < p; j ++) {
    REAL(Rprobs)[j] += (1.0 - REAL(samplemargs)[j])* eta/NC;
    if (REAL(Rprobs)[j] > 1.0) REAL(Rprobs)[j] = 1.0;
//    Rprintf("j = %d sample pip %lf corrected pip = %lf \n", j, REAL(samplemargs)[j], REAL(Rprobs)[j]);
    }
  UNPROTECT(1);
}

/* Used in lm_mcmcbas.c*/
void compute_sampleprobs_modelspace_Bernoulli(SEXP modelspace, SEXP modeldim, SEXP Rsampleprobs, SEXP Rprobs, 
                       int nModels, int p)
{
  int m, j, *model; 
  int *modelVec;
  modelVec = ivecalloc(p);
  memset(modelVec, 0, p * sizeof(int));
  
    for(m=0; m < nModels; m++) {
    memset(modelVec, 0, p * sizeof(int));
    model = INTEGER(VECTOR_ELT(modelspace,m));
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) {
      modelVec[model[j]] = 1.0;
    }
    REAL(Rsampleprobs)[m] = compute_sample_probs_bernoulli(Rprobs, modelVec, p);

  }
}

/* Used in glm_deterministic.c */
void compute_margprobs_old(Bit **models, SEXP Rmodelprobs, double *margprobs, int k, int p)
{
  int m, j;
  double *modelprobs;
  modelprobs = REAL(Rmodelprobs);

  for (j=0; j< p; j++) {
    margprobs[j] = 0.0;
   for(m=0; m< k; m++) {
     if (models[m][j])
        margprobs[j] += modelprobs[m];
      }
    }
}

int no_prior_inclusion_is_1(int p, double *probs) {

  int noInclusionIs1 = 0;
  // loop starts from 1 since the intercept is corrected for in the model prior functions
  for (int i = 1; i < p; i++) { 
  	if (probs[i] > (1.0 - DBL_EPSILON)) {
  		noInclusionIs1++;
  	}
  }
  return noInclusionIs1;
}

/* Is this doing the same as the GetModel_all function? */
void model_to_vec(int *model, int p, SEXP Rmodel) {
  int j;
  memset(model, 0, p * sizeof(int));
  
  for (j = 0; j < LENGTH(Rmodel); j++) {
   model[INTEGER(Rmodel)[j]] = 1;
  }
}

/* Commented in lm_mcmcbas.c*/
double compute_sample_probs_bernoulli(SEXP Rprobs, int *model, int p) {
  int j;
  double pigamma = 1.0;
  for (j = 0; j < p; j++) {
    pigamma *= ((double) model[j])*REAL(Rprobs)[j] + (1.0 - ((double) model[j]))*(1.0 -  REAL(Rprobs)[j]);
  }

  return(pigamma);
}

/*Same name as the newer version (uncomment in the future)
double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior, int noInclusionIs1) {
  const char *family;
  double *hyper_parameters, priorprob = 1.0;


  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  // do not reduce p by the number of predictors that are always included
  // Gitub issue # 87
  if (strcmp(family, "Bernoulli") == 0)
    priorprob = Bernoulli(model, p, hyper_parameters);
  
  // reduce the model space by the number of predictors that are always included 
  p -= noInclusionIs1;
  modeldim -= noInclusionIs1;

  if  (strcmp(family, "Beta-Binomial") == 0)
    priorprob = beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Beta-Binomial") == 0)
    priorprob = trunc_beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Poisson") == 0)
    priorprob = trunc_poisson(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Power-Prior") == 0)
    priorprob = trunc_power_prior(modeldim, p, hyper_parameters);
// Need to add
//  if (strcmp(family, "Hereditary") == 0)
//    priorprob = Hereditary(model, p, hyper_parameters);
  return(priorprob);
}
*/

/*Same name as the newer version (uncomment in the future)
double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior, int noInclusionIs1) {
  const char *family;
  double *hyper_parameters, priorprob = 1.0;


  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  // do not reduce p by the number of predictors that are always included
  // Gitub issue # 87
  if (strcmp(family, "Bernoulli") == 0)
    priorprob = Bernoulli(model, p, hyper_parameters);
  
  // reduce the model space by the number of predictors that are always included 
  p -= noInclusionIs1;
  modeldim -= noInclusionIs1;

  if  (strcmp(family, "Beta-Binomial") == 0)
    priorprob = beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Beta-Binomial") == 0)
    priorprob = trunc_beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Poisson") == 0)
    priorprob = trunc_poisson(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Power-Prior") == 0)
    priorprob = trunc_power_prior(modeldim, p, hyper_parameters);
// Need to add
//  if (strcmp(family, "Hereditary") == 0)
//    priorprob = Hereditary(model, p, hyper_parameters);
  return(priorprob);
}
*/

/* I´m using this version, as changes need to be made in bas_glmFC.R 
   This issue will be dealt with in the future                       */
double compute_prior_probs(int *model, int modeldim, int p, SEXP modelprior) {
  const char *family;
  double *hyper_parameters, priorprob = 1.0;


  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));

  // do not reduce p by the number of predictors that are always included
  // Gitub issue # 87
  if (strcmp(family, "Bernoulli") == 0)
    priorprob = Bernoulli(model, p, hyper_parameters);
  
  // reduce the model space by the number of predictors that are always included 
  // p -= noInclusionIs1;
  // modeldim -= noInclusionIs1;

  if  (strcmp(family, "Beta-Binomial") == 0)
    priorprob = beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Beta-Binomial") == 0)
    priorprob = trunc_beta_binomial(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Poisson") == 0)
    priorprob = trunc_poisson(modeldim, p, hyper_parameters);
  if  (strcmp(family, "Trunc-Power-Prior") == 0)
    priorprob = trunc_power_prior(modeldim, p, hyper_parameters);
// Need to add
//  if (strcmp(family, "Hereditary") == 0)
//    priorprob = Hereditary(model, p, hyper_parameters);
  return(priorprob);
}

/* When hyper = 0.5 for all variables, it´s the familiar uniform model prior*/
double Bernoulli(gsl_vector *index, //int *model,
                 int p, double *hyper) {
  double prior; /* Prior Model Probability*/
  int j;

  /* Loop skips the intercept */
  for (j=1, prior=1.; j < p; j++) {

    int model_j = (int) gsl_vector_get(index, j); // Binary inclusion status
    
    /* hyper [j] = PriorIP for variable j */
    switch(model[j]) {
      case 0:
        prior *= (1. - hyper[j]);
        break;
      case 1:
        prior *= hyper[j]; 
        break;
    }
  }
  return(prior);
}

/* Log-scale evaluation 1st for stability */

/* modeldim = d + 1; hyper = {a, b} */
double beta_binomial(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtact 1 from each */
  
  return(exp(lbeta((double) modeldim - 1.0 + hyper[0], (double) (p - modeldim) + hyper[1]) -
	     lbeta(hyper[0], hyper[1])));
}

double trunc_beta_binomial(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtact 1 from each */

  double prior;
  if ((double) (modeldim -1) <= hyper[2]) {
      prior = exp(lbeta((double) modeldim - 1.0 + hyper[0], (double) (p - modeldim) + hyper[1]) -
		  lbeta(hyper[0], hyper[1]));
      //      Rprintf("pass \n");

    }
  else {prior = 0.0;}

  /*  Rprintf("prior %lf pmodel= %d s0 = %lf\n", prior, modeldim,
      hyper[2]); */
  return(prior);
}

double trunc_poisson(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtract 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
      prior = exp(dpois(modeldim - 1, hyper[0], 1) - ppois(hyper[1], hyper[0], 1, 1) - lchoose((double) p-1, (double) modeldim - 1)); 
    }
  else {prior = 0.0;}

  return(prior);
}

double trunc_power_prior(int modeldim, int p, double *hyper) {
  /* modeldim and p include the intercept so subtract 1 from each */

  double prior = 0.0;
  if ((double) (modeldim -1) <= hyper[1]) {
    prior = exp(-((double) modeldim - 1.0)*((double) hyper[0])*log((double) hyper[1]+1) -
                  (log1mexp((double)(hyper[1]+1)*(double) hyper[0]*log((double) (hyper[1]+1))) - log1mexp((double) hyper[0]*log( (double) (hyper[1]+1)))) -
                  lchoose((double) p-1, (double) modeldim - 1));
    }
  else {prior = 0.0;}

  return(prior);
}

// Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
// This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
// License text is available at https://www.gnu.org/licenses/gpl-3.0.html
// SPDX-License-Identifier: GPL-3.0
//
#include "bas.h"

// Computes posterior model probabilities for each model
void compute_modelprobs (SEXP Rmodelprobs,  SEXP Rlogmarg, SEXP Rpriorprobs, int k)
{ // k stands for the number of models
	int m;
	double nc, bestmarg, *modelprobs, *logmarg, *priorprobs;

  // They are passed as R objects (SEXP), but they are converted to C arrays
	logmarg = REAL(Rlogmarg); //  Input Vector that contains the log marginal likelihoods
	modelprobs = REAL(Rmodelprobs); //  Output Vector that will store the posterior model probabilities
	priorprobs = REAL(Rpriorprobs); // Input Vector that contains the prior model probabilities

	bestmarg = logmarg[0]; //  Initializing the best marginal likelihood to the first model
	nc = 0.0; //  Initializing the normalization constant

	for (m = 0; m < k; m++) { // Finds the maximum log marginal likelihood
		if (logmarg[m] > bestmarg) bestmarg = logmarg[m];
	}

	for (m = 0; m < k; m++) {
    modelprobs[m] = logmarg[m] - bestmarg; /* log (marg / max(marg)) 
    it is basically the renormalized log marginal likelihoods */
		nc += exp(modelprobs[m])*priorprobs[m]; // Denominator in Bayes Theorem for posterior model probabilities
	}

  /* When the entire model space is enumerated, this is the same thing as the original Bayes Theorem.
     When it´s not, this was the way they found the make the normalizing constant comparable among different 
     number of visited / sampled models. */

	for (m = 0; m < k; m++) {
    // From Bayes Theorem, we have that:
	  /*		modelprobs[m] = exp(modelprobs[m]) +
			log(priorprobs[m]) - log(nc)); */
    /* David´s explanation (direct application of the theorem):
       modelprobs[m] = exp(modelprobs[m])) * priorprobs[m] / nc; 
    */
		modelprobs[m] = exp(modelprobs[m] - log(nc))*priorprobs[m]; // Applying Bayes Theorem for posterior model probabilities computation
	}
}

// Computes marginal inclusion probabilities for each predictor (for factors, this means at the level of the levels)
// Input:
//   modelspace   - a list (length k) of integer vectors (each vector = one model's variable indices) 
//   modeldim     - integer vector (length k) giving the number of variables in each model (includes the intercept)
//   Rmodelprobs  - numeric vector (length k) of posterior probabilities for each model
//   margprobs    - (output) C array of length p, filled with marginal inclusion probabilities
//   k            - number of models
//   p            - number of total predictors
void compute_margprobs (SEXP modelspace, SEXP modeldim, SEXP Rmodelprobs, 
                        double *margprobs, int k, int p)
{
	int m, j, *model;
	double *modelprobs;
  // Convert the R object of model posterior probabilities to a C pointer
	modelprobs = REAL(Rmodelprobs);
  // Initialize the (output) C array of length p with marginal inclusion probabilities to 0
	for (j=0; j < p; j++)  margprobs[j] = 0.0;

  // Loop through all models
	for(m=0; m< k; m++) {
		model = INTEGER(VECTOR_ELT(modelspace, m)); // Get the model's variable indices
		// modeldim[m] tells how many variables are in model m (including the intercept)
    for (j = 0; j < INTEGER(modeldim)[m]; j ++) { // Loop through all variables in the model (model is a list of indices)
			// Adding the Posterior Probabilities of this particular model to the PIP for the variable
      margprobs[model[j]] += modelprobs[m]; // model[j] is index of the j-th variable in the model (this index will indeed correspond to a variable)
		}
	}

}


// Used in glm_deterministic.c
// models is a (C 2D binary) matrix of 0/1 indicators of which variables are active in each model instead of a list of indices 
// Arguments:
//   models       - A 2D binary array (Bit**) of size [k x p]
//                  where models[m][j] == 1 if variable j is in model m and p contains the intercept
void compute_margprobs_old(Bit **models, SEXP Rmodelprobs, double *margprobs, int k, int p)
{
  int m, j;
  double *modelprobs;
  // Convert R vector to C array of posterior probabilities
  modelprobs = REAL(Rmodelprobs);

  // Loop through all variables
  for (j=0; j< p; j++) {
    margprobs[j] = 0.0; // Initialize the PIP for the variable to 0
   for(m=0; m< k; m++) { // Loop through all models
     if (models[m][j]) // If the variable is active in the model
        margprobs[j] += modelprobs[m]; // Adding the Posterior Probabilities of this particular model to the PIP for the variable
      }
    }
}

// Returns the number of predictors (excluding the intercept) that have a prior inclusion probability (used for sampling) equal to 1
// Essentially, this counts the number of forced variables (variables that are always included in the model), because we set their initial probability to 1
int no_prior_inclusion_is_1(int p, double *probs) {

  int noInclusionIs1 = 0; // Initialize the counter to 0
  // loop starts from 1 since the intercept is corrected for in the model prior functions
  for (int i = 1; i < p; i++) { 
  	if (probs[i] > (1.0 - DBL_EPSILON)) { // If the prior inclusion probability is close to 1
  		noInclusionIs1++; // Increment the counter
  	}
  }
  return noInclusionIs1; // Returns the final count
}


/* Added / Modified by David */
/*----------------------------------------------------------
 *  compute_prior_probs_MCMC()
 *
 *  Returns the prior model probability for a given model
 *  (represented by the int* 'model' bitmask).
 *
 *  Arguments
 *  ---------
 *    model      : int*  – binary vector of length p (1 = var in model) (includes the intercept)
 *    modeldim   : int   – number of 1's in 'model'
 *    p          : int   – total # of candidate predictors (includes the intercept...)
 *    modelprior : SEXP  – R list with elements
 *                         $family            (string)
 *                         $hyper.parameters  (numeric vector)
 *---------------------------------------------------------*/
double compute_prior_probs_MCMC (int * model, int p, SEXP modelprior,
                                 gsl_matrix *positions, int nofvars,
                                 gsl_vector *levels, gsl_matrix *costs, int n_obs)
 {
  
  const char *family;
  double *hyper_parameters, priorprob = 1.0;

  /*--------------------------------------------------*/
  /*  Pull the prior family name and its hyper-params */
  /*--------------------------------------------------*/
  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));
  

  /* Mudar as prior probs functions para que não seja 
     preciso mudar de um int para um gsl_vector    */

  // Allocate a GSL vector of size p
	gsl_vector *index = gsl_vector_alloc (p);
	// Copy values from model to index
	for (size_t i = 0; i < p; ++i) {
		gsl_vector_set (index, i, (double)model[i]);
	}
  
  gsl_vector_view index_sub = gsl_vector_subvector (index, 1, p - 1);
  gsl_vector *index_ptr = &index_sub.vector;

  if  (strcmp(family, "SBSB") == 0)
    priorprob = SBSB_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "SBC") == 0)
    priorprob = SBC_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "CC") == 0)
    priorprob = CC_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "Uniform") == 0)
    priorprob = C_prior (nofvars, levels);
  if  (strcmp(family, "SB") == 0)
    priorprob = SB_prior (index_ptr, positions, nofvars, levels, p - 1);

  // Marginal Costs for each variable (diagonal of the costs [, -0] matrix)
  double *marginal_costs = (double *) malloc(nofvars * sizeof(double));
  // Initializing the baseline cost
  double c0 = 1000;
  double c0_discount = 1000;         /* high value, just for testing */

  int num_vars_group0 = 0, group_id;
  for (int i = 0; i < nofvars; i++) {
    
    group_id = (int) round (gsl_matrix_get (costs, i, 0)); // 1st Column
    if (group_id == 0)
      num_vars_group0++;

    marginal_costs[i] = gsl_matrix_get (costs, i, i + 1); 
    if (marginal_costs[i] < c0) {
      c0 = marginal_costs[i]; // Updating the baseline cost
    }

    for (int j = 1; j <= nofvars; ++j) {     // j ≤ nofvars, not < 
        double v = gsl_matrix_get(costs, i, j);
        if (v < c0_discount) c0_discount = v;
    }
    
  }

  /* ACHO QUE AGORA JÁ PODEREI ELIMINAR ISTO... */
  if  (strcmp(family, "FND") == 0)
    priorprob = FND_prior (index_ptr, positions, nofvars, levels, p - 1, hyper_parameters, marginal_costs, c0, n_obs);

  if (num_vars_group0 != nofvars) // i.e., less than nofvars...
    c0 = c0_discount;

  if  (strcmp(family, "FNDConst") == 0)
    priorprob = FNDConst_prior (index_ptr, positions, nofvars, levels, p - 1, hyper_parameters, costs, marginal_costs, c0, n_obs);
  if  (strcmp(family, "FNDSB") == 0)
    priorprob = FNDSB_prior (index_ptr, positions, nofvars, levels, p - 1, hyper_parameters, costs, marginal_costs, c0, n_obs);
  
  free (marginal_costs);

  return(priorprob);

}

/*---------------------------------------------------------------------------
 *  SBSBpriorprob(): same as priorSBSB1 in R (wrongly defined)
 *
 *  Computes the **double Scott-Berger prior** πSBSB(M) for a candidate model M.
 *  ─────────────────────────────────────────────────────────────────────────
 *  • 1st Scott–Berger layer (SB1): assigns equal weight to every *variable*
 *    configuration of size (m1 + m2).  m1 = # numeric vars in M,
 *    m2 = # factors in M.
 *
 *   marginal(\gamma, \tau)  =  (1 / choose(nofvars , m1 + m2)) * 1/(nofvars+1)
 *
 *  • 2nd Scott–Berger layer (SB2): for each *factor* present, assigns equal
 *    weight to every subset of its dummy levels that has size v[i].
 *
 *        πSB2(M)  ∝  ∏factor_i  1 / choose(levels[i] , v[i])
 *
 *  prior:  conditional(\delta | \gamma, \tau) × marginal(\gamma, \tau)
 *                   (everything is computed on the log-scale for stability,
 *                   then exponentiated and reciprocated once at the end).
 *
 *  Arguments
 *  ---------
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positionsx : (length nofvars) binary vector –> 1 if column j is *numeric*
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *---------------------------------------------------------------------------*/
double SBSB_prior (gsl_vector * index,
                     gsl_matrix * positions,
                     int          nofvars,
                     gsl_vector * levels,
                     int          p)
{
    double activeVars   = 0.0;   /* # variables (numeric + factors)    */
    double logDen = log (nofvars + 1.0);

    /*------------------------------------------------------*/
    /*  2) For every parent variable i …                    */
    /*     – compute suma  (# active dummies)               */
    /*     – increment activeVars if the variable appears at all   */
    /*------------------------------------------------------*/
    for (int i = 0; i < nofvars; i++) { // for each variable i (row of positions matrix)

        double suma = 0.0;   /* running tally of active dummies for a given variable i */
        int    j    = 0;
        double Li = gsl_vector_get(levels, i); /* # Levels of the variable */

        /* Walk across the design columns j = 0 … p-1 */
        while (suma < Li && j < p) {
        // the number of active dummies for variable i can´t be higher then its number of levels

            /* Add 1 when column j belongs to variable i *and* is in M (is active) */
            suma += gsl_matrix_get(positions, i, j) *
                    gsl_vector_get(index,  j);

            j++;
        }

        // Conditional layer (inverse proportional to the number of models of that DIMENSION, instead of being of that RANK)
        /* If at least one level (only one option for numeric vars) is active,
         * then this parent (competing) variable is in the model. */
        if (suma > 0.0) { // suma is the number of active dummies for variable i
            activeVars++;
            logDen += log (Li) + gsl_sf_lnchoose (Li, suma);   /* the log of the product is the sum´s log*/
        }
    } // we´re adding 0s for numeric variables, so they don´t contribute to the logDen

    logDen += gsl_sf_lnchoose (nofvars, activeVars); // End of marginal layer

    double priorProb = exp(-logDen);

    return priorProb;
    
}

/*==============================================================================
 *  ConstConstpriorprob() <-> equivalent to priorConstConst1 in R (wrongly defined)
 *  Both the marginal and conditional layers are inversely proportional to the number of models
 * 

 *  Computes the **double-constant (CC) prior** probability πCC(M) for a
 *  candidate model M.  The idea is:
 *
 *      • First layer (variables / marginal): every *parent* variable i
 *        (numeric or factor) is selected with probability ½.
 *
 *          ⇒  π₁(M) =  2^-nofvars                         (constant)
 *
 *      • Second layer (factor levels / conditional):  once a factor is in the model,
 *        any NON-EMPTY subset of its dummy levels is equiprobable.
 *        A factor with Lᵢ levels has (2^{Lᵢ} – 1) possible subsets.
 *
 *          ⇒  π₂(M) =  ∏_{factor i ∈ M} (2^{Lᵢ} – 1)^{-1}
 *
 *      • Numeric variables have only one "level", so their contribution at
 *        the second layer is (2¹–1)=1 and drops out.
 *
 *      • The function works on the log–scale for stability and returns
 *        exp(–log-denominator) so that the value is already the reciprocal
 *        of the accumulated terms.
 *
 *  Arguments
 *  ---------
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positionsx : (length nofvars) binary vector –> 1 if column j is *numeric*
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *==============================================================================*/
double CC_prior (gsl_vector * index,
                            gsl_matrix * positions,
                            int          nofvars,
                            gsl_vector * levels,
                            int          p)
{
    double logDen = log (pow (2.0, nofvars));

    /*=== 1. For every parent variable compute v[i] and m1plusm2 = m1+m2 ====*/
    for (int i = 0; i < nofvars; i++) {

        double suma = 0.0;   /* running count of active dummies for var i */
        int    j    = 0;
        double Li = gsl_vector_get(levels, i); /* # levels   */

        /* Walk across all design columns until either:
         *   – every level of var i accounted for OR
         *   – end of column list                                             */
        while (suma < Li && j < p) {

            /* Add 1 if column j belongs to var i AND is active in model M */
            suma += gsl_matrix_get(positions, i, j) *
                    gsl_vector_get(index,  j);

            j++;
        }

        if (suma > 0) {   /* variable i is in model */
          /* Add  log(2^{Li} – 1)  to the denominator      */
          logDen += log(pow(2.0, Li) - 1.0); // we´re including more models than we should (we should exclude the repeated models)
        } // we´re adding 0s for numeric variables, so they don´t contribute to the logDen

    }

    /* Prior probability is the reciprocal of the accumulated factors     */
    double priorProb = exp(-logDen);

    return priorProb;
}

/*==============================================================================
 *  SBConstpriorprob <-> equivalent to priorSBConst1 in R (wrongly defined)
 *
 *  Computes the prior probability of a model under a **Scott–Berger × Constant**
 *  (SB-Const) prior:
 *
 *    - First (/marginal) layer (parent variable inclusion):     Scott–Berger prior
 *          π₁(M) = (1 / choose(nofvars, m1 + m2)) * * 1/(nofvars+1)
 *
 *    - Second (/conditional) layer (factor level subsets):         Constant prior
 *          π₂(M) = ∏_{factors ∈ M} 1 / (2^Lᵢ - 1)
 *
 *  Numeric variables are only affected by the marginal layer.
 *
 *  Parameters:
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *============================================================================*/
double SBC_prior (gsl_vector * index,
                  gsl_matrix * positions,
                  int          nofvars,
                  gsl_vector * levels,
                  int          p)
{
    double logDen = log (nofvars + 1.0);
    double activeVars = 0.0;   // Total count of parent variables in the model (numeric + factor)

    // Loop through each parent variable i (numeric or factor)
    for (int i = 0; i < nofvars; i++) {
        
        double suma = 0.0;  // Running total of selected dummies for variable i
        int j = 0;
        double Li = gsl_vector_get(levels, i); /* # Levels of the variable */

        // Go through each column j until all levels of variable i are processed
        while ((suma < Li) && (j < p)) {
            // Add 1 if col j belongs to var i AND is selected in the model
            suma += gsl_matrix_get(positions, i, j) * gsl_vector_get(index, j);

            j++;
        }

        // If at least one dummy column of var i is selected, it's in the model (active)
        // Conditional Part
        if (suma > 0.0) {
            activeVars++;
            logDen += log(pow(2.0, Li) - 1.0); // We´re subtracting 1 for the EMPTY model (the variable is not in the model)
        } // Numeric variables contribute 0 to this layer, so it´s okay not to ignore them

    }

    logDen += gsl_sf_lnchoose (nofvars, activeVars);

    /* Prior probability is the reciprocal of the accumulated factors     */
    double priorProb = exp(-logDen);

    return priorProb;
}

/* The standard constant prior, inversely proportional to the number of models */
double C_prior (int          nofvars,
                gsl_vector * levels)
{
    double logDen = 0.0;

    /*=== 1. For every competing variable ===*/
    for (int i = 0; i < nofvars; i++) {

        double Li = gsl_vector_get(levels, i); /* # levels */
        logDen += Li * log(2.0);

    }

    /* Prior probability is the reciprocal of the accumulated factors     */
    double priorProb = exp(-logDen);

    return priorProb;
}

/* The standard SB prior where prob over models is inversely 
   proportional to the number of models of that dimension */
double SB_prior (gsl_vector * index,
                 gsl_matrix * positions,
                 int          nofvars,
                 gsl_vector * levels,
                 int          p)
{

    double model_dim = 0.0; // Model dimension (number of active variables / dummies)
    double sum_Li = 0.0; // Sum of the levels for all variables
    double logDen = 0.0;

    for (int i = 0; i < nofvars; i++) { // for each variable i (row of positions matrix)

        double suma = 0.0; /* running tally of active dummies for a given variable i */
        int    j    = 0;
        double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
        
        sum_Li += Li; // Sum of the levels for all variables

        /* Walk across the (p) design columns j = 0 … p-1 */
        while (suma < Li && j < p) {
        // the number of active dummies for variable i can´t be higher then its number of levels

            /* Add 1 when column j belongs to variable i *and* is in M (is active) */
            suma += gsl_matrix_get (positions, i, j) *
                    gsl_vector_get (index,  j);

            j++;
        }

        model_dim += suma;
    
    }

    logDen += log (sum_Li + 1.0) + gsl_sf_lnchoose (sum_Li, model_dim); 
    return exp(-logDen);
      
}

double penalty_cost_func (double cost, double c0, double * b, char *cost_type) {
  
  /* Important to dereference the pointer b (otherwise we´re passing the address of the pointer)
     to get the value of b, *b. */

  if (strcmp (cost_type, "ECP") == 0) {
    return ((pow((cost/c0), *b) - 1));
  }
  else if (strcmp (cost_type, "LCP") == 0) {
    return (((*b * (cost-c0) + c0) / c0));
  }
  else {
    printf("Error: Invalid cost type\n");
    return 0.0;
  }

}

double FND_prior (gsl_vector * index, gsl_matrix * positions, int nofvars, 
                  gsl_vector * levels, int p, double * b, 
                  double * marginal_costs, double c0, int n)
{

  double logDen = 0.0, logNum = 0.0;

  for (int i = 0; i < nofvars; i++) {
      
      double Li = gsl_vector_get (levels, i); /* # Levels of the variable */
      double suma = 0.0; /* running tally of active dummies for a given variable i */
      int j = 0;
      
      /* Page 5 of the Flexible Cost-Penalized Bayesian Model Selection paper*/
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP"); 
      logDen += Li * log (1 + pow (n, -0.5 * cost_penalty));
      
      /* Walk across the (p) design-matrix columns j = 0 … p-1 */
      while (suma < Li && j < p) {
      // the number of active dummies for variable i can´t be higher then its number of levels

          /* Add 1 when column j belongs to variable i *and* is in M (is active) */
          suma += gsl_matrix_get (positions, i, j) *
                  gsl_vector_get (index,  j);

          j++;
      }

      /* If suma is 0, we´re not penalizing the prior probability (variable is not in the model)*/
      /* For categorical variables, incorporated all*/ 
      logNum += suma * cost_penalty;
      
  }

  logNum *= log (n) * (-0.5); 

  double priorProb = exp (logNum - logDen);
  return priorProb;

}

double FNDConst_prior (gsl_vector * index, gsl_matrix * positions, int nofvars, gsl_vector * levels,
                       int p, double * b, gsl_matrix * costs, double * marginal_costs, double c0, int n)
{
  
  // Variable Groups
  int *which_group = malloc(nofvars * sizeof(int));
  
  int max_groupID = 0, num_no_group_variables = 0;
  double logDen_conditional = 0.0, logDen_noGroup = 0.0, logNum_noGroup = 0.0;
  
  for (int i = 0; i < nofvars; i++) {
    
    which_group [i] = (int) round(gsl_matrix_get (costs, i, 0)); // 1st Column
    if (which_group [i] > max_groupID) {
      max_groupID = which_group [i];
    }

    double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
    double suma = 0.0; /* running tally of active dummies for a given variable i */
    int j = 0;
          
    /* Walk across the (p) design-matrix columns j = 0 … p-1 */
    while (suma < Li && j < p) {
    // the number of active dummies for variable i can´t be higher then its number of levels

        /* Add 1 when column j belongs to variable i *and* is in M (is active) */
        suma += gsl_matrix_get (positions, i, j) *
                gsl_vector_get (index,  j);

        j++;
    }
    
    /* INDEPENDENCE SETTING */
    if (which_group [i] == 0) {
      
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP"); 
      logDen_noGroup += log (1 + pow (n, -0.5 * cost_penalty));
    
      // If at least one dummy column of var i is selected, it's in the model (active)
      if (suma > 0.0) { 
        logNum_noGroup += log (n) * cost_penalty * (-0.5); // Costs are accounted for at the marginal level only  
      } 

      num_no_group_variables++;

    }
    
    /* CONDITIONAL PART (see SBConst prior for more info) */
    if (suma > 0) {   /* variable i is in model */
      /* Add  log(2^{Li} – 1)  to the denominator      */
      logDen_conditional += log(pow(2.0, Li) - 1.0); // we´re including more models than we should (we should exclude the repeated models)
    } // we´re adding 0s for numeric variables, so they don´t contribute to the logDen

  }

  double priorProb = exp ((logNum_noGroup - logDen_noGroup) - logDen_conditional);
  
  if (num_no_group_variables == nofvars) {

    free (which_group);
    return priorProb;

  }

  // Initializing with 0
  int *vars_per_group = calloc (max_groupID + 1, sizeof *vars_per_group);

  for (int i = 0; i < nofvars; ++i) {
      int g = which_group[i];
      vars_per_group[g]++;
  }
   
  double sum_logPrior_groups = 0.0; 
  for (int g = 1; g <= max_groupID; g++) { /* Skipping Group 0 */
    
    if (vars_per_group[g] == 0)
        continue;

    size_t m = vars_per_group[g];

    gsl_matrix *pos_buf   = gsl_matrix_alloc (m, p);
    gsl_vector *lev_buf   = gsl_vector_alloc (m);
    gsl_matrix *cost_buf  = gsl_matrix_alloc (m, m);

    size_t row_idx = 0;
    for (int i = 0; i < nofvars; ++i)
        if (which_group[i] == g) {

            /* copy its position row (p columns) */
            for (int j = 0; j < p; ++j) 
                gsl_matrix_set(pos_buf, row_idx, j, gsl_matrix_get(positions, i, j));
            
            /* copy its level */
            gsl_vector_set (lev_buf, row_idx, gsl_vector_get(levels, i));

            size_t col_idx = 0;
            for (int j = 0; j < nofvars; ++j)
                if (which_group[j] == g) {
                    gsl_matrix_set(cost_buf, row_idx, col_idx,
                                   gsl_matrix_get(costs, i, j + 1));
                    col_idx++;
                }

            row_idx++;
        }

    /* It´s in log terms */
    sum_logPrior_groups += prior_group (index, pos_buf, (int) m, lev_buf,
                                        p, b, cost_buf, c0, n);  
    
    gsl_matrix_free (pos_buf); gsl_vector_free (lev_buf); gsl_matrix_free (cost_buf);

  }

  free (which_group); free (vars_per_group);
  
  priorProb *= exp (sum_logPrior_groups);  
  return priorProb;

}

double FNDSB_prior (gsl_vector * index, gsl_matrix * positions, int nofvars, gsl_vector * levels,
                    int p, double * b, gsl_matrix * costs, double * marginal_costs, double c0, int n)
{
  
  // Variable Groups
  int *which_group = malloc(nofvars * sizeof(int));
  
  int max_groupID = 0, num_no_group_variables = 0;
  double logDen_conditional = 0.0, logDen_noGroup = 0.0, logNum_noGroup = 0.0;
  
  for (int i = 0; i < nofvars; i++) {
    
    which_group [i] = (int) round(gsl_matrix_get (costs, i, 0)); // 1st Column
    if (which_group [i] > max_groupID) {
      max_groupID = which_group [i];
    }

    double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
    double suma = 0.0; /* running tally of active dummies for a given variable i */
    int j = 0;
          
    /* Walk across the (p) design-matrix columns j = 0 … p-1 */
    while (suma < Li && j < p) {
    // the number of active dummies for variable i can´t be higher then its number of levels

        /* Add 1 when column j belongs to variable i *and* is in M (is active) */
        suma += gsl_matrix_get (positions, i, j) *
                gsl_vector_get (index,  j);

        j++;
    
    }
    
    /* INDEPENDENCE SETTING */
    if (which_group [i] == 0) {
      
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP");  
      logDen_noGroup += log (1 + pow (n, -0.5 * cost_penalty));
    
      // If at least one dummy column of var i is selected, it's in the model (active)
      if (suma > 0.0) { 
        logNum_noGroup += cost_penalty; // Costs are accounted for at the marginal level only  
      } 

      num_no_group_variables++;

    }
    
    /* CONDITIONAL PART / layer: (inverse proportional to the number 
    of models of that DIMENSION, instead of being of that RANK) */
       
    if (suma > 0.0) {
      // Numeric variables contribute 0 to this layer, so it´s okay not to ignore them
      logDen_conditional += log (Li) + gsl_sf_lnchoose (Li, suma); // the log of the product is the sum´s log
    }

  }
  
  logNum_noGroup *= log (n) * (-0.5);

  double priorProb = exp ((logNum_noGroup - logDen_noGroup) - logDen_conditional);
  
  if (num_no_group_variables == nofvars) {
   
    free (which_group);
    return priorProb;

  }

  // Initializing with 0
  int *vars_per_group = calloc (max_groupID + 1, sizeof *vars_per_group);

  for (int i = 0; i < nofvars; ++i) {
      int g = which_group[i];
      vars_per_group[g]++;
  }
   
  double sum_logPrior_groups = 0.0; 
  for (int g = 1; g <= max_groupID; g++) { /* Skipping Group 0 */
    
    if (vars_per_group[g] == 0)
        continue;

    size_t m = vars_per_group[g];

    gsl_matrix *pos_buf   = gsl_matrix_alloc (m, p);
    gsl_vector *lev_buf   = gsl_vector_alloc (m);
    gsl_matrix *cost_buf  = gsl_matrix_alloc (m, m);

    size_t row_idx = 0;
    for (int i = 0; i < nofvars; ++i)
        if (which_group[i] == g) {

            /* copy its position row (p columns) */
            for (int j = 0; j < p; ++j) 
                gsl_matrix_set(pos_buf, row_idx, j, gsl_matrix_get(positions, i, j));
            
            /* copy its level */
            gsl_vector_set (lev_buf, row_idx, gsl_vector_get(levels, i));

            size_t col_idx = 0;
            for (int j = 0; j < nofvars; ++j)
                if (which_group[j] == g) {
                    gsl_matrix_set(cost_buf, row_idx, col_idx,
                                   gsl_matrix_get(costs, i, j + 1));
                    col_idx++;
                }

            row_idx++;
        }

    /* It´s in log terms */
    sum_logPrior_groups += prior_group (index, pos_buf, (int) m, lev_buf,
                                        p, b, cost_buf, c0, n);  
    
    gsl_matrix_free (pos_buf); gsl_vector_free (lev_buf); gsl_matrix_free (cost_buf);

  }

  free (which_group); free (vars_per_group);
  
  priorProb *= exp (sum_logPrior_groups);  
  return priorProb;

}

// Won´t work if we have variables already forced in to be in the model
double compute_prior_probs_enumeration  (gsl_vector *index, int p, SEXP modelprior,
                                         gsl_matrix *positions, int nofvars,
                                         gsl_vector *levels, gsl_matrix *costs, int n_obs)
{
  
  const char *family;
  double *hyper_parameters, priorprob = 1.0;

  /*--------------------------------------------------*/
  /*  Pull the prior family name and its hyper-params */
  /*--------------------------------------------------*/
  family = CHAR(STRING_ELT(getListElement(modelprior, "family"),0));
  hyper_parameters = REAL(getListElement(modelprior,"hyper.parameters"));
 

  gsl_vector_view index_sub = gsl_vector_subvector (index, 1, p - 1);
  gsl_vector *index_ptr = &index_sub.vector;

  // (alterar eventualmente os argumentos destas funções)
  if  (strcmp(family, "SBSB") == 0)
    priorprob = SBSB_enum_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "SBC") == 0)
    priorprob = SBC_enum_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "CC") == 0)
    priorprob = CC_enum_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "Uniform") == 0)
    priorprob = Constant_prior (index_ptr, positions, nofvars, levels, p - 1);
  if  (strcmp(family, "SB") == 0)
    priorprob = ScottBerger_prior (index_ptr, positions, nofvars, levels, p - 1);

  // Marginal Costs for each variable (diagonal of the costs [, -0] matrix)
  double *marginal_costs = (double *) malloc(nofvars * sizeof(double));
  // Initializing the baseline cost
  double c0_discount = 1000;         /* high value, just for testing */
  double c0 = 1000;

  int num_vars_group0 = 0, group_id;
  for (int i = 0; i < nofvars; i++) {
    
    group_id = (int) round (gsl_matrix_get (costs, i, 0)); // 1st Column
    if (group_id == 0)
      num_vars_group0++;

    marginal_costs[i] = gsl_matrix_get (costs, i, i + 1); 
    if (marginal_costs[i] < c0) {
      c0 = marginal_costs[i]; // Updating the baseline cost
    }

    for (int j = 1; j <= nofvars; ++j) {     // j ≤ nofvars, not < 
        double v = gsl_matrix_get(costs, i, j);
        if (v < c0_discount) c0_discount = v;
    }
    
  }
  
  // FND_prior, FND_enum_prior 
  /* FND_enum não é uma distribuição de probabilidades*/
  if  (strcmp(family, "FND") == 0)
    priorprob = FND_prior (index_ptr, positions, nofvars, levels, p - 1, 
                           hyper_parameters, marginal_costs, c0, n_obs);

  if (num_vars_group0 != nofvars) // i.e., less than nofvars...
    c0 = c0_discount;

  if  (strcmp(family, "FNDConst") == 0)
    priorprob = FNDConst_enum_prior (index_ptr, positions, nofvars, levels, p - 1, hyper_parameters, 
                                     costs, marginal_costs, c0, n_obs);
  if  (strcmp(family, "FNDSB") == 0)
    priorprob = FNDSB_enum_prior (index_ptr, positions, nofvars, levels, p - 1, hyper_parameters, 
                                  costs, marginal_costs, c0, n_obs);

  free (marginal_costs);

  return (priorprob);
}

/*==============================================================================
 *  Constant_prior <-> equivalent to priorConst2 in R
 *
 *  Parameters:
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *============================================================================*/
double Constant_prior (gsl_vector * index,
                       gsl_matrix * positions,
                       int          nofvars,
                       gsl_vector * levels,
                       int          p)      
{
    double logDen = 0.0;

    /*=== 1. For every competing variable */
    for (int i = 0; i < nofvars; i++) {

        double suma = 0.0;   /* running count of active dummies for var i */
        int    j    = 0;
        double Li = gsl_vector_get(levels, i); // # Levels of the variable

        /* Walk across all design columns until either:
         *   – every level of var i accounted for OR
         *   – end of column list                                             */
        while (suma < Li && j < p) {

            /* Add 1 if column j belongs to var i AND is active in model M */
            suma += gsl_matrix_get(positions, i, j) *
                    gsl_vector_get(index,  j);

            j++;
        }
        
        /* if ((suma == Li - 1)  && (Li > 1.0)) {
          return 0.0;
        } */

        /* Building the log-denominator of the Constantprior */
        /* The log of the product is the sum of the logs */
        if (Li > 1.0) {  // Categorical variables
          logDen += log(pow(2.0, Li) - Li);
        }
        else { // Numeric variables
          logDen += log(2.0);
        }
        
    }

    /* Prior probability is the reciprocal of the accumulated factors     */
    double priorProb = exp(-logDen);

    return priorProb;
}

/*==============================================================================
 *  CC_enum_prior <-> equivalent to priorConstConst2 in R
 *
 *  Parameters:
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *============================================================================*/
double CC_enum_prior  (gsl_vector * index,
                       gsl_matrix * positions,
                       int          nofvars,
                       gsl_vector * levels,
                       int          p)      
{

  double logDen  = log (pow (2.0, nofvars)); 
  
  /*=== 1. For every competing variable */
  for (int i = 0; i < nofvars; i++) {

      double suma = 0.0;   /* running count of active dummies for var i */
      int    j    = 0;
      double Li = gsl_vector_get(levels, i); // # Levels of the variable

      /* Walk across all design columns until either:
        *   – every level of var i accounted for OR
        *   – end of column list                                             */
      while (suma < Li && j < p) {

          /* Add 1 if column j belongs to var i AND is active in model M */
          suma += gsl_matrix_get(positions, i, j) *
                  gsl_vector_get(index,  j);

          j++;
      }

      // If the variable is categorical (and is active), we adjust logDen accordingly
      if (Li > 1.0 && suma > 0) {
        if (suma == Li - 1) { // If the factor is saturated, return 0
          return 0.0;
        }    
        logDen += log (pow(2.0, Li) - Li - 1);
        
      }
  } 
  
  /* Prior probability is the reciprocal of the accumulated factors     */
  double priorProb = exp(-logDen);

  return priorProb;        
        
}

/*==============================================================================
 *  SBC_enum_prior <-> equivalent to priorSBConst2 in R
 *
 *  Parameters:
 *    p          : # columns in the overparameterized design matrix (intercept already removed)
 *    index      : (length p) binary vector –> 1 if column j is in M (binary expression of the model)
 *    positions  : (nofvars × p) 0/1 matrix –> row i flags which columns belong
 *                  to variable i (single 1 for numeric, many 1's for a factor)
 *    nofvars    : total # "parent" (competing) variables (numeric + factors)
 *    levels     : (length nofvars) – # levels in each variable (1 for numeric, > 1 for factors)
 *============================================================================*/
double SBC_enum_prior (gsl_vector * index,
                       gsl_matrix * positions,
                       int          nofvars,
                       gsl_vector * levels,
                       int          p)      
{
 
  double activeVars = 0.0;   // Total count of parent variables in the model (numeric + factor)
  double logDen  = log (nofvars + 1.0);  

  /*=== 1. For every competing variable */
  for (int i = 0; i < nofvars; i++) {

      double suma = 0.0;   /* running count of active dummies for var i */
      int    j    = 0;
      double Li = gsl_vector_get (levels, i); /* # Levels of the variable */

      /* Walk across all design columns until either:
        *   – every level of var i accounted for OR
        *   – end of column list                                             */
      while (suma < Li && j < p) {

          /* Add 1 if column j belongs to var i AND is active in model M */
          suma += gsl_matrix_get(positions, i, j) *
                  gsl_vector_get(index,  j);

          j++;
      }

      if (suma > 0) { // If the variable is active, we increment activeVars
        if (suma == Li - 1) { // If the factor is saturated, return 0 
          return 0.0;  /* suma here is always >= 1 so this doesn´t happen for numeric variables */
        }
        if (Li > 1.0) { // If the variable is categorical (and is active), we adjust logDen accordingly
          logDen += log (pow(2.0, Li) - Li - 1);
        }
        activeVars++;
      }

    }
  
  logDen += gsl_sf_lnchoose (nofvars, activeVars);

  /* Prior probability is the reciprocal of the accumulated factors     */
  double priorProb = exp(-logDen);

  return priorProb;        

} 
      
/*==============================================================================
 *  my_choose <-> equivalent to my.choose in R
 *
 *  Auxiliary function to compute the number of models of a given rank 
 *  (see Definition S.1 from professor´s paper Supplementary Materials)
 * 
 *  Parameters:
 *    l          : number of levels of a given factor
 *    j          : number of active levels of a given factor
 *============================================================================*/
double my_choose (int l,
                  int j) 
{

    if (j > l || j < 0 || l <= 0) {
        error ("Bad arguments on my.choose function");
    }

    // Core logic
    if (j < (l - 1)) {
        return choose (l, j);  // Using R's choose()
    } else {
        return 1.0;
    }
}

/*---------------------------------------------------------------------------
 * rank_levels  –>  C equivalent of rank.levels (R)
 * 
 * To be used together with SBSB_enum_prior
 * 
 * Returns a vector with values equal to the number of models of each rank,
 * where indeed the given rank serves as the index of the vector, so, 
 * m2 <= rank <= (sum_i (l_i) | suma > 0, Li > 1) - m2 <=> 
 * 0 <= rank - m2 <= (sum_i (l_i) | suma > 0, Li > 1) - 2*m2, to start indexing from 0.
 * We assume that there´s no saturated model (that situation is accounted for
 * in the SBSB_enum_prior function)
 *
 * Arguments:
 *    L  : array with Li of the *active* factors (length m2)  
 *    m2 : number of active factors in the model
 *
 * Returns: 
 *      An array  result[0 .. len-1]  where
 *
 *          len = Σ_j L[j] - 2*m2 + 1
 *
 *      result[q]  holds  the number of (non-saturated) models whose total
 *      number of active levels is
 *
 *          r = m2 + q    (  m2 ≤ r ≤ Σ_j L_j − m2  )
 *
 *      The caller owns the memory and must  free()  the returned pointer.
 *
 *  Notes
 *  -----
 *  • The last entry (q = len-1) is always set to 1, because there is
 *    exactly one over-parameterised model for  r = Σ_j L_j − m2.
 *---------------------------------------------------------------------------*/
double *rank_levels (double * L,
                     int      m2)
{

    /* ---------- total number of (a₁,…,a_m2) combinations ---------------- */
    long long tot_rows = 1; /* may be large, so 64-bit (p.levels in R, nrow of mm matrix) */
    int sumL = 0;

    for (int j = 0; j < m2; ++j) {
        tot_rows *= (L[j] - 1);      /* Π_j (l_j − 1)                      */
        sumL     +=  L[j];  /* L_j1 + ... + L_jm2 */
    }

    /* admissible r range                                                   */
    // m2 is the minimum admissible rank (n.levels in R)
    const int r_max = sumL - m2;
    const int len   = r_max - m2 + 1; // number of admissible ranks (check formula 15 from the professor´s paper)

    double *result = (double *)calloc(len, sizeof(double)); // Allocating memory for the result vector and initializing it to 0
    if (!result) return NULL; // If the allocation failed, return NULL
    // We free the result vector in SBSB_enum_prior, after calling this (rank_levels) function

    /* ───────────────────────── enumerate combinations ─────────────────────────
    * Allocate a workspace vector `a` that will hold the current tuple
    * (a₁, …, a_{m2}) while we iterate through EVERY combination of
    * active levels.
    *
    *   • `m2`            : number of factors being varied here
    *   • sizeof(int)     : bytes per int on this platform
    *   • total bytes     : m2 * sizeof(int)
    *
    * Each a_j will take integer values in 1, … , (l_j − 1) because we erase
    * the saturated models from the analysis, so we can treat the overparameterized
    * situations as if we were in a saturated one, given a given parametrization.
    * ───────────────────────────────────────────────────────────────────────── */
    int *a = (int *)malloc(m2 * sizeof(int)); /* current (a₁,…,a_{m2}) tuple (data structure)
    for each possible combination of ranks for active factors, this will store the current ranks */

    for (long long row = 0; row < tot_rows; ++row) {

        /* decoding each row into the mixed-radix digits (Lj − 1) -> each factor has its base
        we do this to enumerate the cartesian product 1:(l₁ − 1) × … × 1:(l_m2 − 1), i.e, the set
        of all ordered combinations, one could say */
        /* works also for the special situation of one two-level factor (m2 = 1 and L[0] = 2)*/
        long long tmp = row;
        int r  = 0;                  /* Σ_j a_j                             */
        double cp = 1.0;             /* Π_j my.choose(l_j, a_j)                */

        for (int j = m2 - 1; j >= 0; --j) {
             
            int Lj = (int) round(L[j]);
            int base = Lj - 1;
            a[j] = (int)(tmp % base) + 1;   /* 1 … L_j−1 */
            // Recall that % is the modulo operator (remainder of an integer division)
            tmp /= base;
            
            r  += a[j]; // model rank for this combination
            cp *= my_choose (Lj, a[j]); // round() is used to avoid floating point errors (converting a double to an int)
        }

        result[r - m2] += cp;     /* accumulate for this rank */
    }

    free(a);

   /* --------------------------------------------------------------------------
    CORRECT THE LAST entry, r_max = len - 1, to be 1.0 (the overparameterized model)
   ---------------------------------------------------------------------------*/
    result [len - 1] = 1.0; /* From the repeated models, we keep only the overparameterized... */
    return result;
              
}

/*---------------------------------------------------------------------------
 *  SBSB_enum_prior  –  C equivalent of priorSBSB2 (R)
 *
 *  Arguments (same meaning as before)
 *    index      : p length binary (double) vector, with values 1 if column j is in the model, 0 otherwise
 *    positions  : nofvars × p binary matrix, with values 1 if column j belongs to variable i
 *    nofvars    : (numeric + factor parents)
 *    levels     : 1 × nofvars, with values Li (=1 numeric, ≥2 factor)
 *    p          : # design columns in the overparameterized design matrix (no intercept)
 *---------------------------------------------------------------------------*/
double SBSB_enum_prior (gsl_vector *index,
                        gsl_matrix *positions,
                        int         nofvars,   
                        gsl_vector *levels,
                        int         p)
{
    int    activeVars = 0;
    int    nFactors   = 0;   /* how many Li > 1 (for both active and unactive variables) */
    double sum_levels_act_factors = 0.0;
    int model_rank = 0;

    /* First count the factors so we can allocate small scratch arrays       */
    for (int i = 0; i < nofvars; ++i)
        if (gsl_vector_get(levels, i) > 1.0) ++nFactors;

    double *levelsf      = calloc(nFactors, sizeof *levelsf);     /* will store Li for active factors only */
    int     m2           = 0;                                     /* factor index (can´t be i because nFactors is not the same as nofvars) 
                                                                    + will store the number of active factors in the model */
  
    for (int i = 0; i < nofvars; ++i) {

        double suma = 0.0;
        int j = 0;
        double Li = gsl_vector_get(levels, i);      /* 1 (num) or ≥2  (factor)  */

        /* Walk across the design columns j = 0 … p-1 */
        while (suma < Li && j < p) {
        // the number of active dummies for variable i can´t be higher then its number of levels

            /* Add 1 when column j belongs to variable i *and* is in M (is active) */
            suma += gsl_matrix_get(positions, i, j) *
                    gsl_vector_get(index,  j);

            j++;
        }

        if (suma > 0.0) {
          
          // If Li = 1 (numeric), it will never be equal to suma + 1, as suma can´t be 0 (making this only "work" for factors)
          if (suma == Li - 1) { // If the factor is saturated, return 0 (this model doesn´t exist)
              free (levelsf); 
              return 0.0;
            }

          activeVars++; /* parent is in M */

          if (Li > 1.0) {

            levelsf [m2]     = Li; // Stores the number of levels for active factors
            ++m2; // Stores the number of active factors in the model

            model_rank += suma; // Stores the model rank (not considering numeric variables)
            sum_levels_act_factors += Li; // l_j1 + ... + l_jm2 (regardless of whether the model is saturated or not)

            if (suma == Li) {
              model_rank--; // we´re removing 1 for oversaturated factors
            } 
          }
        }                                                
    }

    double logDen = log(nofvars + 1.0) + gsl_sf_lnchoose(nofvars, activeVars); /* Marginal layer */

    /* ---------- *no* factor present (constant marginal layer) only ---- */
    if (m2 == 0) {
        free (levelsf);
        return exp(-logDen); /* Prior Model Probability */
    }

    /* ----------  FACTORS ARE PRESENT  -------------------------------------- */
    // Obtain the vector with the number of models for each admissible rank:
    double *num_models_same_rank = rank_levels (levelsf, m2); // m2 serves to "control the dimension" of the levelsf vector
    /* Given a vector of the levels for the active factors (l1,l2,...,l_m2) and its length, this function 
    computes how many models (not saturated) there are with the same number of active levels (r) 
    such that m2 <= r <= (sum_i (l_i) | suma > 0) - m2, i.e. 0 <= r - m2 <= (sum_i (l_i) | suma > 0) - 2 * m2 */
    free (levelsf);

    /* index in the 0-based C array is r - m2, where r is the model rank */
    int r_minus_m2 = model_rank - m2;

    /* -------------------------------------------------------------------------------------
     Just to check if everything is okay
    ------------------------------------------------------------------------------------- */
    /* max_value_for_r_minus_m2 = Σ l_j – 2·m2 */
    double max_r_minusm2_value = sum_levels_act_factors - 2.0 * (double)m2;
    /* Safety check: r_minus_m2 (index) must lie inside the vector we just received */
    if (r_minus_m2 < 0 || r_minus_m2 > (int)max_r_minusm2_value) {              
        free (num_models_same_rank);
        error ("index (r - m2) out of bounds");   
    }

    /* Adding the Conditional Layer (Formula 15 from the professor´s paper) */
    logDen += log (sum_levels_act_factors - 2.0 * m2 + 1.0) + 
              log (num_models_same_rank[r_minus_m2]); 

    free (num_models_same_rank);

    double priorProb = exp(-logDen);
    return priorProb;

}	

double *rank_levels2 (double * arg,
                      int nofvars, 
                      int len_out) {

  /* ---------- total number of (a₁,…,a_nofvars) combinations ---------------- */
    long long tot_rows = 1; /* may be large, so 64-bit (p.levels in R, nrow of mm matrix) */
    
    for (int j = 0; j < nofvars; ++j) {
        tot_rows *= arg[j];      /* Π_j ((l_j - 1) - 0 + 1) = Π_j (l_j) ; more rows than in rank_levels */
        /* We are multiplying by 2 for both numeric variables and factors with two levels,
           meaning that we can either include them in a model or not */
        /* We did a trick in the ScottBerger_prior to attribute two levels to numeric variables
           in order for this to work */
    }

    double *result = (double *)calloc(len_out, sizeof(double)); // Allocating memory for the result vector and initializing it to 0
    if (!result) return NULL; // If the allocation failed, return NULL
    // We free the result vector in ScottBerger_prior, after calling this (rank_levels) function

    /* ───────────────────────── enumerate combinations ─────────────────────────
    * Allocate a workspace vector `a` that will hold the current tuple
    * (a₁, …, a_{nofvars}) while we iterate through EVERY combination of
    * active levels.
    *
    *   • `nofvars`       : number of competing variables being varied here
    *   • sizeof(int)     : bytes per int on this platform
    *   • total bytes     : nofvars * sizeof(int)
    *
    * Each a_j will take integer values in 0, … , (l_j − 1) because we erase
    * the saturated models from the analysis, so we can treat the overparameterized
    * situations as if we were in a saturated one, given a given parametrization.
    * ───────────────────────────────────────────────────────────────────────── */
    int *a = (int *)malloc(nofvars * sizeof(int)); /* current (a₁,…,a_{nofvars}) tuple (data structure)
    for each possible combination of ranks for active factors, this will store the current ranks */

    for (long long row = 0; row < tot_rows; ++row) {

        /* decoding each row into the mixed-radix digits (arg_j) -> each variable has its base
        we do this to enumerate the cartesian product 0:(l₁ − 1) × … × 0:(l_{nofvars} − 1), i.e,
        the set of all ordered combinations, one could say */
        // NEED TO CONFIRM THIS!!!
        /* works also for the special situation of one two-level factor (nofvars = 1 and arg[0] = 2)
        => i can i have 2 models (1 of rank 1 and another of rank 0)*/
        /* Also works for the null model */

        long long tmp = row;
        int r  = 0;                  /* Σ_j a[j] = Combination´s model rank     */
        double cp = 1.0;             /* Π_j my.choose(arg_j, a [j])                */

        for (int j = nofvars - 1; j >= 0; --j) {
            
            int arg_j = (int) round(arg[j]); // converting a double to an int
            // arg_j will act as the base for the mixed-radix digits
            // int base = arg_j;
            a[j] = (int)(tmp % arg_j);   /* 0 … (arg_j − 1) */
            // Recall that % is the modulo operator (remainder of an integer division)
            tmp /= arg_j;
            
            r  += a[j]; // model rank for this combination
            cp *= my_choose (arg_j, a[j]); 
        }

        result[r] += cp;     /* accumulate for this rank */
    }

    free(a);

   /* --------------------------------------------------------------------------
    CORRECT THE LAST entry, max_rank = len_out - 1, to be 1.0 (the overparameterized model)
   ---------------------------------------------------------------------------*/
    result [len_out - 1] = 1.0; /* From the repeated models, we keep only the overparameterized... */
    /* This is the situation where all competing factors are saturated / oversaturated 
    and all the numeric variables are included in the model */
    return result;
              
}

/*---------------------------------------------------------------------------
 *  ScottBerger_prior  –  C equivalent of priorSB2 (R)
 *
 *  Under this function, prior probabalities are inversely proportional to the number of models
 *  of a given rank (and constant over all dimensions).
 * 
 *  For repeated models, we keep only the overparameterized one
 *
 *  Arguments (same meaning as before)
 *    index      : p length binary (double) vector, with values 1 if column j is in the model, 0 otherwise
 *    positions  : nofvars × p binary matrix, with values 1 if column j belongs to variable i
 *    nofvars    : (numeric + factor parents)
 *    levels     : 1 × nofvars, with values Li (=1 numeric, ≥2 factor)
 *    p          : # design columns in the overparameterized design matrix (no intercept)
 *---------------------------------------------------------------------------*/
double ScottBerger_prior (gsl_vector * index, gsl_matrix * positions, int nofvars, gsl_vector * levels, int p)
{

  int model_rank = 0;
  int max_rank = 0;

  /* Creating an array of nofvars doubles on the heap, initiali<ing all of them to 0,
   and keeping its address in the pointer arg. */
  double * arg = calloc (nofvars, sizeof *arg); 
 
  /*=== 1. For every competing variable (row in the positions matrix)*/
  for (int i = 0; i < nofvars; i++) {

    double suma = 0.0;   /* running count of active dummies for var i */
    int    j    = 0;
    double Li = gsl_vector_get (levels, i); // Levels of the variable (>1 for factors, 1 for numeric)
    arg [i] = fmax (Li, 2.0); // Trick to make the rank_levels2 function also work for numeric variables
    max_rank += arg [i] - 1; // Sums 1 for both numeric variables and factors with two levels

    /* Walk across all design columns until either:
      *   – every level of var i accounted for OR
      *   – end of column list                                             */
    while (suma < Li && j < p) {

        /* Add 1 if column j belongs to var i AND is active in model M */
        suma += gsl_matrix_get (positions, i, j) *
                gsl_vector_get (index,  j);

        j++;
    }
    
    if (suma > 0.0) {
          
      // If Li = 1 (numeric), it will never be equal to suma + 1, as suma can´t be 0 (making this only "work" for factors)
      if (suma == Li - 1) { /* Erase from the set of competing models the saturated models */
          free (arg); 
          return 0.0;
        }


      if (Li > 1.0) {

        model_rank += round(suma); // Stores the model rank (not considering numeric variables)
        // round() is used to avoid floating point errors (converting a double to an int)

        if (suma == Li) {
          model_rank--; // we´re removing 1 for oversaturated factors (not full-rank)
        } 

      } 
      else {
        model_rank += 1; // suma can´t exceed 1 for numeric variables
      }
    }

  }

  double logDen = log (max_rank + 1);
  /* Essentially, it´s the number of different admissible ranks a model can have:
   * min_rank = 0 <= r <= max_rank = sum (L_i - 1) + k  
     => Length of admissible ranks = max_rank - min_rank + 1 = max_rank + 1 */

   /* -------------------------------------------------------------------------------------
  The Scott and Berger prior is constant over all dimensions and inversely proportional to 
  number of models of that given rank (there´s a single model of rank 0 and of maximum /
  full rank).
  ------------------------------------------------------------------------------------- */
  
  /* ff contains the matrix with the count of models of each rank for each admissible rank,
     with 0 <= r <= sum (L_i - 1) + k */
  double *ff = rank_levels2 (arg, nofvars, max_rank + 1);
  free (arg);
  
  /* Adjusting with the number of models of a given admissible rank... */
  logDen += log(ff[model_rank]); 
  free (ff);

  /* Prior probability is the reciprocal of the accumulated factors     */
  double priorProb = exp(-logDen);
  return priorProb;

}

double FND_enum_prior (gsl_vector *index, gsl_matrix *positions, int nofvars, gsl_vector *levels,
                       int p, double *b, double *marginal_costs, double c0, int n) 
{

  double logDen = 0.0, logNum = 0.0;

  for (int i = 0; i < nofvars; i++) {
      
      double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
      double suma = 0.0; /* running tally of active dummies for a given variable i */
      int j = 0;
      
      /* Page 5 of the Flexible Cost-Penalized Bayesian Model Selection paper*/
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP"); 
      logDen += Li * log (1 + pow (n, -0.5 * cost_penalty));
      
      /* Walk across the (p) design-matrix columns j = 0 … p-1 */
      while (suma < Li && j < p) {
      // the number of active dummies for variable i can´t be higher then its number of levels

          /* Add 1 when column j belongs to variable i *and* is in M (is active) */
          suma += gsl_matrix_get (positions, i, j) *
                  gsl_vector_get (index,  j);

          j++;
      }
      
      /* It´s as if we´re creating restrictions to the model space... */
      /* PROBLEM: The prior probabilities won´t add up to 1 over the entire model space */
      /* Check comments before priorFND2 in R */
      if (suma > 0.0 && suma == Li - 1) {
        return 0.0; // Saturated models
      }

      /* If suma is 0, we´re not penalizing the prior probability (variable is not in the model)*/
      /* For categorical variables, penalize the dummies included */
      logNum += suma * cost_penalty;
      
  }

  logNum *= log (n) * (-0.5);

  double priorProb = exp (logNum - logDen);
  return priorProb;

}

/* Here costs is a square matrix (nofvars x nofvars) */
double prior_group (gsl_vector * index, gsl_matrix * positions, int nofvars, 
                    gsl_vector * levels, int p, double * b, gsl_matrix * costs,
                    double c0, int n)
{
  
  if (nofvars == 0) {
    return 0.0; /* No prior probability assigned to this group */
    /* Leave as log 1 = 0.0, so there´s no problem doing the product of prior probabilities */
  }

  /* Reduces to the independence setting when nofvars is equal to 1 */
  // Marginal Costs for each variable (diagonal of the costs [, -0] matrix)
  int *active_vars = (int *) malloc(nofvars * sizeof(int));

  for (int s = 0; s < nofvars; s++) { /* I´ve corrected nofvars in the call */
      
      double Ls = gsl_vector_get(levels, s); /* # Levels of the variable */
      double suma = 0.0; /* running tally of active dummies for a given variable i */
      int r = 0;
            
      /* Walk across the (p) design-matrix columns r = 0 … p-1 */
      while (suma < Ls && r < p) {
      // the number of active dummies for variable i can´t be higher then its number of levels

          /* Add 1 when column j belongs to variable i *and* is in M (is active) */
          suma += gsl_matrix_get (positions, s, r) *
                  gsl_vector_get (index,  r);

          r++;
      }

      active_vars [s] = (int) suma;

  }

  double priorprob = 0.0;
  for (int s = 0; s < nofvars; s++) { // nofvars is now the number of variables of this group

    /* Page 5 of the Flexible Cost-Penalized Bayesian Model Selection paper*/
    /* Using the exponential cost prior (ECP) */
    double cost_penalty = penalty_cost_func (gsl_matrix_get (costs, s, s), c0, b, "ECP");
    double aux_logDen = log (1 + pow (n, -0.5 * cost_penalty)); double aux_logNum = 0.0;

    if (active_vars [s] > 0) {
      aux_logNum = log (n) * cost_penalty * (-0.5); // Costs are accounted for at the marginal level only 
    }
    
    for (int j = 0; j < nofvars; j++) {
        
        if (j == s) {
          continue;
        }
        
        double custo = c0; /* Just to check; initialize with 0*/
        if (active_vars [s] > 0) {
          custo = gsl_matrix_get (costs, j, s);
        }
        else {
          custo = gsl_matrix_get (costs, j, j);
        }
        
        double new_penalty = penalty_cost_func (custo, c0, b, "ECP");

        aux_logDen += log (1 + pow (n, -0.5 * new_penalty));
      
        if (active_vars [j] > 0) {
          aux_logNum += log (n) * new_penalty * (-0.5);
        }

    }

    double add = exp (aux_logNum - (aux_logDen + log (nofvars)));
    priorprob += add;

  }   
      
  free (active_vars);

  return log (priorprob);

}


double FNDConst_enum_prior (gsl_vector *index, gsl_matrix *positions, int nofvars, gsl_vector *levels,
                            int p, double *b, gsl_matrix *costs, double * marginal_costs, double c0, int n) 
{ 
  
  // Variable Groups
  int *which_group = malloc(nofvars * sizeof(int));

  int num_no_group_variables = 0, max_groupID = 0;
  double logDen_conditional = 0.0, logDen_noGroup = 0.0, logNum_noGroup = 0.0;
  
  for (int i = 0; i < nofvars; i++) {
    
    which_group [i] = (int) round(gsl_matrix_get (costs, i, 0)); // 1st Column
    if (which_group [i] > max_groupID) {
      max_groupID = which_group [i];
    }

    double suma = 0.0; /* running tally of active dummies for a given variable i */
    int j = 0;
    double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
          
    /* Walk across the (p) design-matrix columns j = 0 … p-1 */
    while (suma < Li && j < p) {
    // the number of active dummies for variable i can´t be higher then its number of levels

        /* Add 1 when column j belongs to variable i *and* is in M (is active) */
        suma += gsl_matrix_get (positions, i, j) *
                gsl_vector_get (index,  j);

        j++;
    }
    
    /* INDEPENDENCE SETTING */
    if (which_group [i] == 0) {
      
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP"); 
      logDen_noGroup += log (1 + pow (n, -0.5 * cost_penalty));
    
      // If at least one dummy column of var i is selected, it's in the model (active)
      if (suma > 0.0) { 
        logNum_noGroup += log (n) * cost_penalty * (-0.5); // Costs are accounted for at the marginal level only  
      } 

      num_no_group_variables++;

    }

    /* CONDITIONAL PART */
    if (suma > 0) { // If the variable is active, we increment activeVars
      if (suma == Li - 1) { // If the factor is saturated, return 0 
        free (which_group);
        return 0.0;  /* suma here is always >= 1 so this doesn´t happen for numeric variables */
      }
      if (Li > 1.0) { // If the variable is categorical (and is active), we adjust logDen accordingly
        logDen_conditional += log (pow(2.0, Li) - Li - 1);
      }
    }

  }
 
  double priorProb = exp ((logNum_noGroup - logDen_noGroup) - logDen_conditional);
  
  if (num_no_group_variables == nofvars) {

    free (which_group);
    return priorProb;

  }

  // Initializing with 0
  int *vars_per_group = calloc (max_groupID + 1, sizeof *vars_per_group);

  for (int i = 0; i < nofvars; ++i) {
      int g = which_group[i];
      vars_per_group[g]++;
  }
   
  double sum_logPrior_groups = 0.0; 
  for (int g = 1; g <= max_groupID; g++) { /* Skipping Group 0 */
    
    // If the group was not seen
    if (vars_per_group[g] == 0) {
      // printf ("Aconteceu-me isto");
      continue;
    }
        
    size_t m = vars_per_group[g];

    gsl_matrix *pos_buf   = gsl_matrix_alloc (m, p);
    gsl_vector *lev_buf   = gsl_vector_alloc (m);
    gsl_matrix *cost_buf  = gsl_matrix_alloc (m, m);

    size_t row_idx = 0;
    for (int i = 0; i < nofvars; ++i)
        if (which_group[i] == g) {

            /* copy its position row (p columns) */
            for (int j = 0; j < p; ++j) 
                gsl_matrix_set(pos_buf, row_idx, j, gsl_matrix_get(positions, i, j));
            
            /* copy its level */
            gsl_vector_set (lev_buf, row_idx, gsl_vector_get(levels, i));

            size_t col_idx = 0;
            for (int j = 0; j < nofvars; ++j)
                if (which_group[j] == g) {
                    gsl_matrix_set(cost_buf, row_idx, col_idx,
                                   gsl_matrix_get(costs, i, j + 1));
                    col_idx++;
                }

            row_idx++;
        }

    /* It´s in log terms */
    sum_logPrior_groups += prior_group (index, pos_buf, (int) m, lev_buf,
                                        p, b, cost_buf, c0, n);  
    
    gsl_matrix_free (pos_buf); gsl_vector_free (lev_buf); gsl_matrix_free (cost_buf);

  }

  free (which_group); free (vars_per_group);

  priorProb *= exp (sum_logPrior_groups);
  return priorProb;

}

double FNDSB_enum_prior (gsl_vector *index, gsl_matrix *positions, int nofvars, gsl_vector *levels,
                         int p, double *b, gsl_matrix *costs, double * marginal_costs, double c0, int n) 
{ 
  
  // Variable Groups, Model Rank & # of Competing Factors
  int *which_group = malloc(nofvars * sizeof(int));

  int num_no_group_variables = 0, max_groupID = 0, model_rank = 0, nFactors = 0;
  double logDen_conditional = 0.0, logDen_noGroup = 0.0, logNum_noGroup = 0.0, sum_levels_act_factors = 0.0;

  /* First count the factors so we can allocate small scratch arrays       */
    for (int i = 0; i < nofvars; ++i)
        if (gsl_vector_get(levels, i) > 1.0) ++nFactors;
  
  double *levelsf      = calloc (nFactors, sizeof *levelsf);     /* will store Li for active factors only */
  int     m2           = 0;                                     /* factor index (can´t be i because nFactors is not the same as nofvars) 
                                                                   + will store the number of active factors in the model */

  for (int i = 0; i < nofvars; i++) {
    
    which_group [i] = (int) round(gsl_matrix_get (costs, i, 0)); // 1st Column
    if (which_group [i] > max_groupID) {
      max_groupID = which_group [i];
    }

    double suma = 0.0; /* running tally of active dummies for a given variable i */
    int j = 0;
    double Li = gsl_vector_get(levels, i); /* # Levels of the variable */
          
    /* Walk across the (p) design-matrix columns j = 0 … p-1 */
    while (suma < Li && j < p) {
    // the number of active dummies for variable i can´t be higher then its number of levels

        /* Add 1 when column j belongs to variable i *and* is in M (is active) */
        suma += gsl_matrix_get (positions, i, j) *
                gsl_vector_get (index,  j);

        j++;
    }


    /* CONDITIONAL PART */
    if (suma > 0) { // If the variable is active, we increment activeVars
      
      if (suma == Li - 1) { // If the factor is saturated, return 0 
        free (which_group); free (levelsf);
        return 0.0;  /* suma here is always >= 1 so this doesn´t happen for numeric variables */
      }

      if (Li > 1.0) {

            levelsf [m2]  = Li; // Stores the number of levels for active factors
            ++m2; // Stores the number of active factors in the model

            model_rank += suma; // Stores the model rank (not considering numeric variables)
            sum_levels_act_factors += Li; // l_j1 + ... + l_jm2 (regardless of whether the model is saturated or not)

            if (suma == Li) {
              model_rank--; // we´re removing 1 for oversaturated factors
            } 
      }

    }
    
    /* INDEPENDENCE SETTING */
    if (which_group [i] == 0) {
      
      /* Using the exponential cost prior (ECP) */
      double cost_penalty = penalty_cost_func (marginal_costs[i], c0, b, "ECP"); 
      logDen_noGroup += log (1 + pow (n, -0.5 * cost_penalty));
    
      // If at least one dummy column of var i is selected, it's in the model (active)
      if (suma > 0.0) { 
        logNum_noGroup += log (n) * cost_penalty * (-0.5); // Costs are accounted for at the marginal level only  
      } 

      num_no_group_variables++;

    }

  }
  
  double priorProb = exp (logNum_noGroup - logDen_noGroup);


  if (m2 != 0) {
    /* ----------  FACTORS ARE PRESENT  -------------------------------------- */
    // Obtain the vector with the number of models for each admissible rank:
    double *num_models_same_rank = rank_levels (levelsf, m2); // m2 serves to "control the dimension" of the levelsf vector
    /* Given a vector of the levels for the active factors (l1,l2,...,l_m2) and its length, this function 
    computes how many models (not saturated) there are with the same number of active levels (r) 
    such that m2 <= r <= (sum_i (l_i) | suma > 0) - m2, i.e. 0 <= r - m2 <= (sum_i (l_i) | suma > 0) - 2 * m2 */
    free (levelsf);
    
    /* index in the 0-based C array is r - m2, where r is the model rank */
    int r_minus_m2 = model_rank - m2;

    /* max_value_for_r_minus_m2 = Σ l_j – 2·m2 */
    double max_r_minusm2_value = sum_levels_act_factors - 2.0 * (double)m2;
    /* Safety check: r_minus_m2 (index) must lie inside the vector we just received */
    if (r_minus_m2 < 0 || r_minus_m2 > (int)max_r_minusm2_value) {              
        free (num_models_same_rank);
        error ("index (r - m2) out of bounds");   
    }

    /* Adding the Conditional Layer (Formula 15 from the professor´s paper) */
    logDen_conditional += log (sum_levels_act_factors - 2.0 * m2 + 1.0) + 
                          log (num_models_same_rank[r_minus_m2]); 

    free (num_models_same_rank);

    priorProb *= exp (-logDen_conditional);

  }

  
  if (num_no_group_variables == nofvars) {

    free (which_group); 
    return priorProb;

  }

  // Initializing with 0
  int *vars_per_group = calloc (max_groupID + 1, sizeof *vars_per_group);

  for (int i = 0; i < nofvars; ++i) {
      int g = which_group[i];
      vars_per_group[g]++;
  }
   
  double sum_logPrior_groups = 0.0; 
  for (int g = 1; g <= max_groupID; g++) { /* Skipping Group 0 */
    
    // If the group was not seen
    if (vars_per_group[g] == 0) {
      // printf ("Aconteceu-me isto");
      continue;
    }
        
    size_t m = vars_per_group[g];

    gsl_matrix *pos_buf   = gsl_matrix_alloc (m, p);
    gsl_vector *lev_buf   = gsl_vector_alloc (m);
    gsl_matrix *cost_buf  = gsl_matrix_alloc (m, m);

    size_t row_idx = 0;
    for (int i = 0; i < nofvars; ++i)
        if (which_group[i] == g) {

            /* copy its position row (p columns) */
            for (int j = 0; j < p; ++j) 
                gsl_matrix_set(pos_buf, row_idx, j, gsl_matrix_get(positions, i, j));
            
            /* copy its level */
            gsl_vector_set (lev_buf, row_idx, gsl_vector_get(levels, i));

            size_t col_idx = 0;
            for (int j = 0; j < nofvars; ++j)
                if (which_group[j] == g) {
                    gsl_matrix_set(cost_buf, row_idx, col_idx,
                                   gsl_matrix_get(costs, i, j + 1));
                    col_idx++;
                }

            row_idx++;
        }

    /* It´s in log terms */
    sum_logPrior_groups += prior_group (index, pos_buf, (int) m, lev_buf,
                                        p, b, cost_buf, c0, n);  
    
    gsl_matrix_free (pos_buf); gsl_vector_free (lev_buf); gsl_matrix_free (cost_buf);

  }

  free (which_group); free (vars_per_group);

  priorProb *= exp (sum_logPrior_groups);
  return priorProb;

}