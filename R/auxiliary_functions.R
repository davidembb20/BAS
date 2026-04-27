# R equivalent to compute_modelprobs in C 
# Similar to .renormalize.postprobs in bas_glmFactors.R, but this was created by me
compute_posterior <- function (logmarg, priorprobs) 
{

    # To use with MCMC, when renormalize = TRUE
    # It uses just the unique models to compute the posterior model probabilities

    # Maximum log marginal likelihood
	best_logmarg <- max (logmarg) 
    # Normalizing Constant (Bayes´s Theorem denominator)
	nc <- sum (exp (logmarg - best_logmarg) * priorprobs) 
	# Applying Bayes´s Theorem to get the posterior model probabilities
	postprobs <- exp (logmarg - best_logmarg) * priorprobs / nc

	return (postprobs)

}

# R equivalent to SBSB_enum_prior in C (but here for "all" prior probs)
priorSBSB2 <- function (models_matrix, positions) 
{
	
	# The corrected SB-SB prior prob (conditionally: inversely proportional to the number of models of that rank)
	
    levels <- rowSums (positions)
	num_vars <- nrow (positions) # Number of competing variables
	levels_f <- levels [levels > 1] # Stores the levels for factors only (>1)

	new_priorprobs <- rep (0, nrow (models_matrix))
	
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)

        # Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
		# Active Factors (with num of active Levels)
		act.levelsf <- act.vars [levels > 1]

        log_Marg <- log (num_vars + 1) + lchoose (num_vars, sum (act.vars != 0))

		# If the model does not contain any factor:
	    if (sum (act.levelsf) == 0) {
			new_priorprobs [i] <- exp (-log_Marg)
		}
		else {
			# Obtain the vector with the number of models of each rank:
			numberof <- rank.levels (levels_f[act.levelsf > 0]) 
			# The argument for rank.levels is a vector of levels for the active factors only

            # If the model is not "saturated" nor overparameterized
			if (sum(act.levelsf >= (levels_f-1)) == 0){ # => act.levelsf < levelsf - 1, for every factor
				aux <- sum (act.levelsf) # ii) from Result S.2 of the supplementary materials: k(gamma,delta) -> model dimension
				m2 <- length (levels_f[act.levelsf > 0]) # Number of active factors
				log_Cond <- log (numberof[as.character(aux)]) + log (sum (levels_f[act.levelsf > 0]) - 2 * m2 + 1) 
				# Formula (15) from professor´s paper, pp.6
				new_priorprobs [i] <- exp (-(log_Cond + log_Marg))
			}
			# If the model contains at least one factor "saturated" then return 0:
			else if (sum(act.levelsf[act.levelsf>0] == (levels_f[act.levelsf>0]-1)) >= 1) {
				new_priorprobs [i] <- 0
			}
			else { # Overparameterized Models
				aux <- sum (act.levelsf[act.levelsf == levels_f] - 1) + sum(act.levelsf[act.levelsf < levels_f]) # ii) from Result S.2 of the supplementary materials: k(gamma,delta) -> model dimension
				m2 <- length (levels_f[act.levelsf > 0]) # Number of active factors
				log_Cond <- log (numberof[as.character(aux)]) + log (sum (levels_f[act.levelsf > 0]) - 2 * m2 + 1) 
				# Formula (15) from professor´s paper, pp.6
				new_priorprobs [i] <- exp (-(log_Cond + log_Marg))
			}
		}
	}

	return (new_priorprobs)
		
}

# R equivalent to my.choose in C
#  Auxiliary function to compute the number of models of a given rank 
#  (see Definition S.1 from professor´s paper Supplementary Materials)
# 
#  Parameters:
#    L          : number of levels of a given factor
#    j          : number of active levels of a given factor
my.choose <- function (L, j) 
{
	
	if (j > L || j < 0 || L <= 0) {
		stop ("Bad arguments on my.choose\n")
	}
	
	if (j < (L-1)) {
		return (choose (L, j)) 
	} else {
		return (1) 
	}
}

# Expected to be used in combination with SBSB <-> R equivalent to rank_levels in C
rank.levels <- function (levelsf)
{

  # ---------------------------------------------------------------------------
  #  PURPOSE
  # ---------------------------------------------------------------------------
  #  For a set of *active* qualitative factors whose numbers of levels are
  #  given in the vector `levelsf` = (l₁, l₂, …, l_m2) with l_j ≥ 2,
  #  this routine counts how many different (non-saturated) models exist
  #  for every admissible rank, r, of *active levels*:
  #
  #        m2  ≤  r  ≤  Σ_j l_j  − m2
  #        0  ≤ r - m2 ≤  Σ_j l_j  − 2 * m2
  #
  #  A model is described by how many levels a₁, a₂, …, a_p (with
  #  1 ≤ a_j ≤ l_j − 1) are retained (/active) in each factor. The
  #  number of ways to pick exactly a_j levels in factor j is
  #  my.choose (l_j, a_j) (implemented before; check definition S.1.
  #  of the supplementary materials)
  #
  #  The function returns a *named* integer vector: the names are the
  #  admissible r-values, and the entries are the corresponding counts.
  #
  # ---------------------------------------------------------------------------
  #  SPECIAL CASE:  a single two-level factor
  # ---------------------------------------------------------------------------
  if (length (levelsf) == 1 && levelsf [1] == 2) {
    result <- c (1)            # exactly one model (just the overparameterized)
    names (result) <- "1"      # r = 1 rank of the corresponding model
    return (result)
  }
  # ---------------------------------------------------------------------------
  #  ENUMERATE **ALL** COMBINATIONS  (a₁, …, a_m2)  WITH 1 ≤ a_j ≤ l_j − 1
  #  For each factor, the maximum rank is L_j - 1
  # ---------------------------------------------------------------------------
  p.levels <- prod (levelsf - 1)      # |{(a₁,…,a_m2)}| = Π_j (l_j − 1)
  n.levels <- length (levelsf)        # m2 = number of active factors
  
  # Matrix `mm` will hold every combination of active levels of factors row-wise 
  # (i.e., we have a combination for each row) 
  # Actually, the values in each column represent the different possible ranks a factor can have (based on their active levels)
  mm <- matrix (0, nrow = p.levels, ncol = n.levels)
  colnames (mm) <- paste0 ("F", seq_len(n.levels)) # Concatenates the string "F" with the sequence of numbers from 1 to n.levels (m2)
  
  # Fill each column with the appropriate repeating pattern so that
  # the Cartesian product 1:(l₁ − 1) × … × 1:(l_m2 − 1) is produced
  if (n.levels > 1) { # If there are more than one active factor (/column)
    for (i in 1:(n.levels - 1)) { # For each column (factor) except the last one? 
      mm[, i] <- rep (1:(levelsf[i] - 1),
	                  # How many times each value is repeated before moving to the next one
					  # For lower (higher) values of i, the number of repetitions is higher (lower)
                      each = prod(levelsf[(i + 1):n.levels] - 1),
                      length.out = p.levels)
    }
	# Filling the last column
    i <- i + 1
    mm[, i] <- rep(1:(levelsf[i] - 1), length.out = p.levels)
  } else {
	# n.levels = 1 and levelsf[1] != 2 
	# p.levels = L_j - 1 > 1, in this case
    mm [, 1] <- rep(1:(levelsf[1] - 1), length.out = p.levels)
  }
  
  # ---------------------------------------------------------------------------
  #  ADD TWO EXTRA COLUMNS
  #    • sum.act.levels  =  r (rank) = Σ_j a_j
  #    • combin.prod     =  Π_j my.choose(l_j, a_j)   (filled later)
  # ---------------------------------------------------------------------------
  mm <- cbind (mm, sum.act.levels = rowSums (mm), combin.prod = 0)
  # ---------------------------------------------------------------------------
  #  COMPUTE  combin.prod  FOR EVERY ROW in mm
  #  Check Result S.2 iii) of the supplementary materials
  #  Each row is just a parcel of the summation
  #  We´re computing the number of models of that given rank for that given 
  #  combination of factor-ranks; but as we can we reach the same full model 
  #  rank by changing factor-ranks among the different factors, we have to sum
  # combin.prod for those rows...
  # ---------------------------------------------------------------------------
  for (row in seq_len(p.levels)) {
    s <- 1
    for (j in seq_len(n.levels)) {
      s <- s * my.choose(levelsf[j], mm [row, j])
    }
    mm [row, "combin.prod"] <- s
  }
  
  # ---------------------------------------------------------------------------
  #  ACCUMULATE COUNTS FOR EACH ADMISSIBLE  r, rank, possible.values
  # ---------------------------------------------------------------------------
  possible.values <- n.levels:(sum(levelsf) - n.levels) # m2 <= r <= sum(l_j) - m2
  result <- setNames(integer(length(possible.values)), possible.values)
  # Initializing the result vector with 0s and attaching names to each element
  #  of the it using the values in possible.values

  for (r in possible.values) {
    rows.with.r <- which(mm[, "sum.act.levels"] == r) # Extracting the rows with the same rank
    result[as.character(r)] <- sum(mm[rows.with.r, "combin.prod"]) # Summing "combin.prod" for those rows
  }
  
  # ---------------------------------------------------------------------------
  #  CORRECT THE LAST  r  (saturated / oversaturated) → exactly *one* model
  # ---------------------------------------------------------------------------
  result[as.character(max(possible.values))] <- 1 # From the repeated models, we keep only the overparameterized...
  
  return(result)
}

# R equivalent to ScottBerger_prior in C (but here for "all" prior probs)
priorSB2 <- function (models_matrix, positions) 
{
	
    # The corrected SB prior prob (inversely proportional to the number of models of that rank)

    levels <- rowSums (positions)
	# num_vars <- nrow (positions) # Number of competing variables
	levels_f <- levels [levels > 1] # Stores the levels for factors only (>1)
	# ff contains the matrix with the count of models of each rank for each combination of levels / variables
    ff <- rank.levels2 (levels)
	logDen <- log (length (ff))

	new_priorprobs <- rep (0, nrow (models_matrix))
	
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)

        # Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]

		# This prior is only based on the rank of the model:
		model_rank <- sum (pmin(act.levelsf, levels_f - 1)) + sum (act.vars[levels == 1])

		# For the null model:
	    if (sum(model_rank==0)==length(model_rank)) { # Can I change this for model_rank == 0?
			new_priorprobs [i] <- exp (-logDen)
		}
		# If the model contains at least one level saturated then return 0:
		else if (sum(act.levelsf[act.levelsf>0] == (levels_f[act.levelsf>0]-1)) >= 1) {	
			new_priorprobs [i] <- 0
		}
		else { # keep the rest (oversaturated models)
			logDen <- logDen + log (ff[as.numeric(names(ff)) == model_rank])
			new_priorprobs [i] <- exp (-logDen)
		}

	}

	return (new_priorprobs)

}

# R equivalent to rank_levels2 in C
rank.levels2 <- function (levelsfull) 
{
	#Expected to be used in combination with SB functions
	#Given a vector of p factors and k variables levelsfull=(l1,l2,...,lp,1,...,1), 
	#(for covariates we have the ones), this function computes how many
	#models there are of all possible ranks 0<=r<=sum(l_i-1)+k
	#(Because of a trick This works even if some of the levelsfull=1 (so numerical vars)
	#but I created the rank.levels3 (see BayesVarSel) exactly for that purpose
	#but it is much more complex and it does the same)
	
	#the resulting table is the same independently of li=1 or li=2
	levelsfull[levelsfull==1]<- 2 # This mean we will treat all the numerical variables as if
	# they were factors with two levels...

	# If we have a single candidate variable (numeric or a factor with two levels)
	if (length(levelsfull)==1){
		if (levelsfull==2){ 
			result<- c(1,1) # Either we include the variable or we don´t, so we have two possible distinct models...
			names(result)<- c("0","1")
			return(result)
			}
	}
	
	p.levels<- prod(levelsfull); # Defined differently in rank.levels function
	n.levels<- length(levelsfull)
	mm<- matrix(0, nrow=p.levels, ncol=n.levels)
	colnames(mm)<- paste("F", seq_along(levelsfull), sep="")
	
	# Cartesian Product 0:(L₁ - 1) × … × 0:(L_p - 1) × (0:1)^k
	# 0:(L₁ - 1) × … × 0:(L_p - 1) × (0:(2-1))^k
	if (n.levels>1){ # If more than one variable...
		for (i in 1:(n.levels-1)){
			mm[,i]<- rep(0:(levelsfull[i]-1), each=prod(levelsfull[(i+1):n.levels]), length.out=p.levels)
		}
			i<- i+1 # Last column
			mm[,i]<- rep(0:(levelsfull[i]-1), each=1, length.out=p.levels)
		}
	else mm[,1]<- rep(0:(levelsfull[1]-1), each=1, length.out=p.levels)
	
	mm<- cbind(mm, rowSums(mm))
	mm<- cbind(mm, 0)
	colnames(mm)[n.levels + (1:2)]<- c("sum.act.levels", "combin.prod")
	s<- 1
	for (i in 1:p.levels){
		for (j in 1:n.levels){
			s<- s*my.choose(levelsfull[j], mm[i,j])
		}
		mm[i, n.levels+2]<- s; s<- 1
	}

	possible.values<- min(mm[,"sum.act.levels"]):max(mm[,"sum.act.levels"])
	result<- rep(0, length(possible.values))
	names(result)<- possible.values
	for (i in possible.values){
		these<- which(mm[,"sum.act.levels"]==i)
		result[as.character(i)]<- sum(mm[these,"combin.prod"])		
	}

	#rectify for the saturated or oversaturated:
	result[as.character(max(possible.values))] <- 1
	return(result)	
}

# R equivalent to Constant_prior in C (but here for "all" prior probs)
priorConst2 <- function (models_matrix, positions) 
{

    levels <- rowSums (positions)
    # Stores the levels for factors only 
	levels_f <- levels [levels > 1]

    # Constant prior over unique models...
	new_priorprobs <- rep (0, nrow (models_matrix))

	# Log of the Number of Unique Models:
	lrui.number <- sum(log(2^levels - levels))
	# We´re adding log (2) for numeric variables...
	
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)

		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]

		# If the model contains at least one level saturated then return 0:
		if (sum(act.levelsf[act.levelsf>0] == (levels_f[act.levelsf>0]-1)) >= 1) {
			new_priorprobs [i] <- 0
		}
		else { # keep the rest (oversaturated or not saturated)
			new_priorprobs [i] <- exp (-lrui.number)
		}
	}

	return (new_priorprobs)

}

# R equivalent to priorConstConst2 in C (but here for "all" prior probs)
priorConstConst2 <- function (models_matrix, positions) 
{

    # The corrected Const-Const prior prob (conditionally: proportional to a constant for unique models)

    levels <- rowSums (positions)
	num_vars <- nrow (positions) # Number of competing variables
	levels_f <- levels [levels > 1] # Stores the levels for factors only 

	new_priorprobs <- rep (0, nrow (models_matrix))
	
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)

        # Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t (positions)
		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]

        log_Marg <- num_vars * log (2) # Prior prob = 2 ^ (-num_vars)

		# If the model does not contain any factor:
	    if (sum (act.levelsf) == 0) {
			new_priorprobs [i] <- exp (-log_Marg)
		}
		else {

			# If the model contains at least one level saturated then return 0:
			if (sum(act.levelsf[act.levelsf>0] == (levels_f[act.levelsf>0]-1)) >= 1) {
				new_priorprobs [i] <- 0
			}
			else { # keep the rest (oversaturated or not saturated)
				log_Cond <- sum (log (2^levels_f[act.levelsf>0] - levels_f[act.levelsf>0] - 1))
		        # For each active factor, I have to remove the null model (as that would mean
				# the factor wouldn´t be active) and also the case where the factor is 
				# saturated (so the number of active levels = number of levels - 1)
				new_priorprobs [i] <- exp (-(log_Cond + log_Marg))
			}
		}
	}

	return (new_priorprobs)

}

# R equivalent to SBC_enum_prior in C (but here for "all" prior probs)
priorSBConst2 <- function (models_matrix, positions) 
{

	# The corrected SB-Const prior prob (conditionally: proportional to a constant for unique models)

    levels <- rowSums (positions)
	num_vars <- nrow (positions) # Number of competing variables
	levels_f <- levels [levels > 1] # Stores the levels for factors only 

	new_priorprobs <- rep (0, nrow (models_matrix))
	
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)

        # Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]

        log_Marg <- log (num_vars + 1) + lchoose (num_vars, sum (act.vars != 0))

		# If the model does not contain any factor:
	    if (sum (act.levelsf) == 0) {
			new_priorprobs [i] <- exp (-log_Marg)
		}
		else {

			# If the model contains at least one level saturated then return 0:
			if (sum(act.levelsf[act.levelsf>0] == (levels_f[act.levelsf>0]-1)) >= 1) {
				new_priorprobs [i] <- 0
			}
			else { # keep the rest (oversaturated or not saturated)
				log_Cond <- sum (log (2^levels_f[act.levelsf>0] - levels_f[act.levelsf > 0] - 1))
				new_priorprobs [i] <- exp (-(log_Cond + log_Marg))
			}
		}
	}

	return (new_priorprobs)

}

# Check matrix.rank.levels function in BayesVarSel (appeared but wasn´t being used in resamplingConst)
# This function also calls the integer.base.b_C function - check BayesVarSel for details)

# R equivalent to penalty_cost_func in C
cost_penalty <- function (costs, b, c0, cost_type) 
{

  if (cost_type == "ECP") {
    return ((costs / c0)^b - 1)
  } 
  else if (cost_type == "LCP") {
    return (((b * (costs - c0) + c0) / c0))
  } 
  else {
    stop ("Error: Invalid cost type")
  }

}


# R equivalent to FND_enum_prior in C (but here for "all" prior probs)

## Acho que não devíamos usar isto!!!
priorFND2 <- function (models_matrix, positions, costs, b, n) 
{

    levels <- rowSums (positions) 
	levels_f <- levels [levels > 1] # Stores the levels for factors only

	new_priorprobs <- rep (0, nrow (models_matrix))

	marginal_costs <- diag (costs) # Diagonal of the costs matrix (after excluding the first column in the call)
    c0 <- min (marginal_costs)

	penalty <- cost_penalty (marginal_costs, b, c0, "ECP")
    logDen <- sum (levels * log (1 + (n ^ (-0.5 * penalty))))

	for (i in seq_len (nrow (models_matrix))) { # For each model (index)
		
		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
        # Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]
				
	    # If the model contains at least one level saturated then return 0:
		if (sum(act.levelsf[act.levelsf > 0] == (levels_f[act.levelsf > 0] - 1)) >= 1) {
			new_priorprobs [i] <- 0
		}

		logNum <- log (n) * sum (act.vars * penalty * (-0.5))
		new_priorprobs [i] <- exp (logNum - logDen)
	
	}

	return (new_priorprobs)

}

prior_no_Group <- function (models_matrix, positions, costs, c0, b, n) 
{
    new_priorprobs <- rep (0, nrow (models_matrix))
	marginal_costs <- diag (costs) # Diagonal of the costs matrix (after excluding the first column in the call)

	# For this function (prior_no_Group2), c0 is given in the call
	penalty <- cost_penalty (marginal_costs, b, c0, "ECP") 
    logDen <- sum (log (1 + (n ^ (-0.5 * penalty))))
    
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)
		
		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t(positions)
        
		logNum <- log (n) * sum (as.numeric (act.vars > 0) * penalty) * (-0.5)
		new_priorprobs [i] <- exp (logNum - logDen)
	
	}

	return (new_priorprobs)

}

# Quando esta função é chamada já só estamos a trabalhar com variáveis de um mesmo grupo
# It reduces to the independence setting when we only one variable in a group...
priorSpecificGroup <- function (models_matrix, positions, costs, c0, b, n) 
{
   
	num_vars <- nrow (positions) # Number of competing variables
	new_priorprobs <- rep (0, nrow (models_matrix))

	marginal_costs <- diag (costs)  
	# THe discounted costs are the ones that are not in the diagonal of the costs matrix...
	penalty <- cost_penalty (marginal_costs, b, c0, "ECP")

	for (i in seq_len (nrow (models_matrix))) { # For each model (index)
		
		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t (positions)

		Den <- 0; Num <- 0; 

		for (s in seq_len (num_vars)) {

            # j <- s

			auxDen <- log (1 +  (n ^ (-0.5 * penalty[s]))) 
            auxNum <- log (n) * penalty [s] * as.numeric (act.vars [s] > 0) * (-0.5)

            for (j in seq_len (num_vars)) { # For each variable in the group
                
				if (j == s) {
					next # Skip
				}
                
				condition <- marginal_costs[j] - as.numeric (act.vars [s] > 0) * (marginal_costs[j] - costs [j, s])
                new_penalty <- cost_penalty (condition, b, c0, "ECP")
				
				auxDen <- auxDen + log (1 +  (n ^ (-0.5 * new_penalty)))
				auxNum <- auxNum + (log (n) * new_penalty * as.numeric (act.vars [j] > 0)) * (-0.5)

			}

            Den <- Den + exp (-auxDen) 
            Num <- Num + exp (-auxNum)

		}

		new_priorprobs [i] <- Num / Den
	
	}

	return (new_priorprobs / num_vars)

}

# Instead of giving every single dummy a prior cost penalization, I´ll just give it if the 
# (categorical) variable is included in the model...
priorFNDConst2 <- function (models_matrix, positions, costs, b, n) 
{

	group <- costs [, 1] # Group for each variable (1st column of the costs matrix)
	num_groups <- length (unique (group))
    num_models <- nrow (models_matrix)

	costs <- costs [, -1] # Excluding the first column (group)
	if (sum(group > 0) > 0) {
		# There is at least one group (i.e., not just zeros)
		c0 <- min (costs) # Use full marginal-cost matrix
	} else {
		c0 <- min (diag(as.matrix(costs))) # Minimum Marginal Cost
	}

	#	if (is.null(num_models) || is.na(num_models))
	#		stop("models_matrix não tem linhas (nrow é NA).")

	df <- matrix (0, nrow = num_models, ncol = num_groups + 1L)
	# An extra column for the conditional prob...

	for (g in seq_len(num_groups)) {
		idx <- group == (g - 1)

		if (!any(idx)) {                     # <-- no rows in this group
			df[, g] <- 0                     # or NA, or simply next
			next                             # skip to the next g
		}

		positions_group <- positions[idx, ]
		costs_group     <- costs [idx, idx]

		if (g == 1) { # => group = 0
			df[, 1] <- log(prior_no_Group(models_matrix,
										  positions_group,
										  costs_group,
										  c0, b, n))
		} else {
			df[, g] <- log(priorSpecificGroup(models_matrix,
											  positions_group,
											  costs_group,
											  c0, b, n))
		}
	}

    
	# Defining the conditional prior probabilities (Constant over Unique Models)
    levels <- rowSums (positions) 
	levels_f <- levels [levels > 1] # Stores the levels for factors only 
    
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)
		
		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t (positions)
		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]

		if (sum (act.levelsf) == 0) {
			df [i, ncol(df)] <- 0 # There´s no conditional (log 1 = 0)
		}
		else {
			# If the model contains at least one level saturated then return -Inf (we´re working with logs):
			if (sum(act.levelsf[act.levelsf > 0] == (levels_f[act.levelsf > 0] - 1)) >= 1) {
				df [i, ncol(df)] <- -Inf
			}
			else { # keep the rest (oversaturated or not saturated)
				# log Denominator for the Conditional
				df [i, ncol(df)] <- -sum (log (2^levels_f[act.levelsf>0] - levels_f[act.levelsf > 0] - 1))
			}
		}
    
	}

	logpriorprobs <- as.vector (rowSums (df)) # Log prior probabilities for each model
	return (exp (logpriorprobs))

}

priorFNDSB2 <- function (models_matrix, positions, costs, b, n) 
{

	group <- costs [, 1] # Group for each variable (1st column of the costs matrix)
	num_groups <- length (unique (group))
    num_models <- nrow (models_matrix)

	costs <- costs [, -1] # Excluding the first column (group)
	if (sum(group > 0) > 0) {
		# There is at least one group (i.e., not just zeros)
		c0 <- min (costs) # Use full marginal-cost matrix
	} else {
		c0 <- min (diag(as.matrix(costs))) # Minimum Marginal Cost
	}

	
	df <- matrix (0, nrow = num_models, ncol = num_groups + 1)
	# An extra column for the conditional prob...

	for (g in seq_len(num_groups)) {
		idx <- group == (g - 1)

		if (!any(idx)) {                     # <-- no rows in this group
			df[, g] <- 0                     # or NA, or simply next
			next                             # skip to the next g
		}

		positions_group <- positions[idx, ]
		costs_group     <- costs [idx, idx]

		if (g == 1) { # => group = 0
			df[, 1] <- log(prior_no_Group(models_matrix,
										  positions_group,
										  costs_group,
										  c0, b, n))
		} else {
			df[, g] <- log(priorSpecificGroup(models_matrix,
											  positions_group,
											  costs_group,
											  c0, b, n))
		}
	}

	# Defining the conditional prior probabilities (Constant over Unique Models)
    levels <- rowSums (positions) 
	levels_f <- levels [levels > 1] # Stores the levels for factors only 
    
	for (i in seq_len (nrow (models_matrix))) { # For each model (index)
		
		# Active Variables (not at the level of the levels)
		act.vars <- models_matrix [i, ] %*% t (positions)

		# Active Levels for Factors
		act.levelsf <- act.vars [levels > 1]
						
		if (sum (act.levelsf) == 0) {
			df [i, ncol(df)] <- 0 # There´s no conditional (log 1 = 0)
		}
		else {
		
			# Obtain the vector with the number of models of each rank:
			numberof <- rank.levels (levels_f[act.levelsf > 0]) 
			# The argument for rank.levels is a vector of levels for the active factors only
		
			# If the model is not saturated nor oversaturated
			if (sum(act.levelsf >= (levels_f-1)) == 0) { # => act.levelsf < levelsf - 1, for every factor
				aux <- sum (act.levelsf)
				m2 <- length (levels_f[act.levelsf > 0]) # Number of active factors
				# Formula (15) from professor´s paper, pp.6
				logDen_Cond <- log (numberof[as.character(aux)]) + log (sum (levels_f[act.levelsf > 0]) - 2 * m2 + 1) 
				df [i, ncol(df)] <- -logDen_Cond
			}
			# If the model contains at least one level saturated then return -Inf (we´re working with logs):
			else if (sum(act.levelsf[act.levelsf > 0] == (levels_f[act.levelsf > 0] - 1)) >= 1) {
				df [i, ncol(df)] <- -Inf
			}
			else { # keep the rest (oversaturated models)
				aux <- sum (act.levelsf[act.levelsf == levels_f] - 1) + sum (act.levelsf[act.levelsf < levels_f])
				m2 <- length (levels_f[act.levelsf > 0])
				# Formula (15) from professor´s paper, pp.6
				logDen_Cond <- log (numberof[as.character(aux)]) + log (sum (levels_f[act.levelsf > 0]) - 2 * m2 + 1) 
				df [i, ncol(df)] <- -logDen_Cond
			}

		}
    
	}
	
	logpriorprobs <- as.vector (rowSums (df)) # Log prior probabilities for each model
	return (exp (logpriorprobs))

}

model_cost <- function (models_active_vars, costs)
{   

	tot.cost <- numeric (nrow(models_active_vars))

	if (sum (costs [, 1]) == 0.00) { # No groups / no cost structures...
		for (i in seq_len (nrow(models_active_vars))) { # For each model (index)
			tot.cost <- as.double(as.vector (models_active_vars %*% diag (costs [, -1])))
		}
		return (tot.cost)
	}

	group <- costs [, 1] # Group for each variable (1st column of the costs matrix)
	num_groups <- length (unique (group))
	group_ids <- sort(unique (group), decreasing = FALSE)
    num_models <- nrow (models_active_vars)

	costs <- costs [, -1] # Excluding the first column (group)

	df <- matrix (0.00, nrow = num_models, ncol = num_groups)


	for (g in group_ids) {
		
		g_idx <- which (group_ids == g)
		idx <- group == g

		if (!any(idx)) {                     # <-- no rows in this group
			df [, g_idx] <- 0.00                 # or NA, or simply next
			next                             # skip to the next g
		}
        
		costs_group <- costs [idx, idx, drop = FALSE]

		if (g == 0) {
			for (i in seq_len (nrow (models_active_vars))) {
				df [i, 1] <- as.double (sum (models_active_vars[i, idx] * diag (costs_group)))
			}
		} else {
			
			for (i in seq_len (nrow (models_active_vars))) { # For each model (index)
				
				if (sum (models_active_vars [i, idx]) <= 1.00) {
                    df [i, g_idx] <- as.double (sum (models_active_vars[i, idx] * diag (costs_group)))
				}
				else { # More than one variable active in a group (activate cost structures)

					active <- models_active_vars [i, idx] > 0.00 
					df [i, g_idx] <- as.double (min (colSums(costs_group [active, active, drop = FALSE])))
				    # We use the minimum, because we´re not in a sequential decision-making framework,
					# so we can choose the cheapest scenario, if it exists
				}
					
			}
			
		}
	}

	tot.cost <- as.double (as.vector (rowSums (df)))
	return (tot.cost)
}



