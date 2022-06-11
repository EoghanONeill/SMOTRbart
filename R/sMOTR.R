#' @import Rcpp
#' @importFrom mvtnorm 'rmvnorm'
#' @importFrom stats 'rgamma' 'runif' 'dnorm' 'sd' 'rnorm' 'pnorm'
#' @importFrom MCMCpack 'rdirichlet'
#' @importFrom truncnorm 'rtruncnorm'
#' @useDynLib SMOTRbart, .registration = TRUE
#' @export
smotr_bart = function(x,
                     y,
                     sparse = TRUE,
                     vars_inter_slope = TRUE,
                     ntrees = 10,
                     node_min_size = 5,
                     alpha = 0.95,
                     beta = 2,
                     nu = 3,
                     lambda = 0.1,
                     sigma2 = 1,
                     nburn = 1000,
                     npost = 1000,
                     nthin = 1,
                     a = 1,
                     b = 1,
                     tau_b = 1,
                     ancestors = FALSE) {

  x = as.data.frame(x)

  # Quantities needed for prediction
  center = apply(x, 2, mean)
  scale = apply(x, 2, sd)
  aux.X = lapply(x, unique) # Checking how many unique values each variable has
  unique.values.X = unlist(lapply(aux.X, length))

  center[which(unique.values.X<=2)] = 0
  scale[which(unique.values.X<=2)] = 1

  # The original X-matrix is stored as X_orig, this will be useful for the covariate sampling
  X_orig = x

  X = as.matrix(cbind(1,scale(x, center, scale))) # standardising the covariates and adding an intercept

  # The standardized X-matrix is stored as X_stand, whis will be useful for the tree modifications and design matrix construction
  X_stand = X

  # Extract control parameters
  node_min_size = node_min_size

  # Extract MCMC details
  TotIter = nburn + npost*nthin # Total of iterations

  #Initialize logarithmic likelihood
  log_lik = 0

  # Storage containers, in this case a container for the beta's, the variance of beta's, bandwidths and log likelihood is added
  store_size = npost
  tree_store = vector('list', store_size)
  beta_store = vector('list', length = store_size)
  tau_store = vector('list', length = store_size)
  sigma2_store = rep(NA, store_size)
  tau_b_store = rep(NA, store_size)
  y_hat_store = matrix(NA, ncol = length(y), nrow = store_size)
  tree_fits_store = matrix(0, ncol = ntrees, nrow = length(y))
  log_lik_store = rep(NA, store_size)

  # Scale the response target variable
  y_mean = mean(y)
  y_sd = sd(y)
  y_scale = (y - y_mean)/y_sd
  n = length(y_scale)

  # Create a list of trees for the initial stump
  curr_trees = create_stump(num_trees = ntrees,
                            y = y_scale,
                            X = X_stand)

  # Initialise the new trees as current one
  new_trees = curr_trees

  # Initialise the predicted values to zero
  predictions = numeric(n)
  pred_mat = matrix(NA,n,ntrees)

  # Set up a progress bar
  pb = utils::txtProgressBar(min = 1, max = TotIter,
                             style = 3, width = 60,
                             title = 'Running rBART...')


  tau_rate = 1
  tau = rep(tau_rate, ntrees) #sample set of bandwidth for each tree from the exponential distribution

  #Initialize two additional objects of type list to store the beta's and bandwidths

  beta_hat = vector('list', length = ntrees)
  tau_vec = vector('list', length = ntrees)

  # Start the MCMC iterations loop
  for (i in 1:TotIter) {

    utils::setTxtProgressBar(pb, i)

    # If at the right place, store everything
    if((i > nburn) & ((i - nburn) %% nthin) == 0) {
      curr = (i - nburn)/nthin
      tree_store[[curr]] = curr_trees
      beta_store[[curr]] = beta_hat
      tau_store[[curr]] = tau_vec
      sigma2_store[curr] = sigma2
      tau_b_store[curr] = tau_b
      y_hat_store[curr,] = predictions
      log_lik_store[curr] = log_lik
    }


    # Start looping through trees
    for (j in 1:ntrees) {

      X = X_stand #For each tree, the design matrix start as the standardized matrix

      current_partial_residuals = y_scale - predictions + tree_fits_store[,j]

      # The number of covariates to sample from and the sample probabilities as based on the original X-matrix
      p = ncol(X_orig)
      s = rep(1/p,p)

      # Propose a new tree via grow/change/prune/swap
      type = sample(c('grow', 'prune', 'change'), 1)

      if(get_nterminal(curr_trees[[j]]) == 1){ type = 'grow' }# If the tree is a stump then grow it
      if(i < max(floor(0.1*nburn), 5)){ type = 'grow' }# Grow for the tree for the first few iterations

      new_trees[[j]] = update_tree(y = y_scale,
                                   X = X_stand, #Note that the standardized matrix is used for the tree construction
                                   type = type,
                                   curr_tree = curr_trees[[j]],
                                   node_min_size = node_min_size,
                                   s = s)


     # Start calculating the design matrix
     # This object determines this most important characteristics of each terminal nodes
     # such as the parent, move to the left, split variable and value
     anc = get_branch(curr_trees[[j]])

    if(is.null(anc)){
      X = matrix(1, nrow = nrow(X_stand), ncol = 1) #If the ancestors object is null, this means the tree is a stump and the design matrix will be one vector
     }
    else{
      #This function creates the phi-matrix based on the standardized data, ancestors object, the phi-function and the bandwidth
      # phi_matrix = t(apply(X_stand,1,phi, anc = anc, tau = tau[[j]]))

      tautemp = tau[[j]]

      phi_matrix = phi_app(X_stand, anc, tautemp)

      #This function construct our design matrix based on the phi-matrix and the covariates that are included in the tree
      X = design_matrix(X_stand,curr_trees[[j]],phi_matrix)
    }

      # Prior for the beta vector
      p = ncol(X)
      V = rep(1/tau_b, p)
      inv_V = 1/V

      # Compute the log of the marginalized likelihood and log of the tree prior for the current tree
       conditional = conditional_tilde(curr_trees[[j]],
                                X,
                                current_partial_residuals,
                                sigma2,
                                V,
                                inv_V,
                                nu,
                                lambda,
                                tau_b,
                                ancestors,
                                ntrees)
       l_old = conditional + get_tree_prior(curr_trees[[j]], alpha, beta)


      #Now the same principle is applied for constructing the design matrix, now using the new tree
      anc_new = get_branch(new_trees[[j]])

      if(is.null(anc_new)){
        X_new = matrix(1, nrow = nrow(X_stand), ncol = 1) #If the ancestors object is null, this means the tree is a stump and the design matrix will be one vector
      }
      else{
        #This function creates the phi-matrix based on the standardized data, ancestors object, the phi-function and the bandwidth
        # phi_matrix_new = t(apply(X_stand,1,phi, anc = anc_new, tau = tau[[j]]))

        tautemp = tau[[j]]

        phi_matrix_new = phi_app(X_stand, anc_new, tautemp)


        #This function construct our design matrix based on the phi-matrix and the covariates that are included in the tree
        X_new = design_matrix(X_stand,new_trees[[j]],phi_matrix_new)
      }


      # New prior for the beta vector
      p_new = ncol(X_new)
      V_new = rep(1/tau_b, p_new)
      inv_V_new = 1/V_new

      # Compute the log of the marginalized likelihood and the log of the tree prior for the new tree
      conditional_new = conditional_tilde(new_trees[[j]],
                                X_new,
                                current_partial_residuals,
                                sigma2,
                                V_new,
                                inv_V_new,
                                nu,
                                lambda,
                                tau_b,
                                ancestors,
                                ntrees)
     l_new = conditional_new + get_tree_prior(new_trees[[j]], alpha, beta)


      # Compute the alpha fir the Metropolis-Hastings step based on the type of tree modification and the likelihoods
      a = alpha_mh(l_new,l_old, curr_trees[[j]],new_trees[[j]], type)


      if(min(1,a) > runif(1)) { # In case the alpha is bigger than a uniformly sampled value between zero and one

        curr_trees[[j]] = new_trees[[j]] # The current tree "becomes" the new tree, if the latter is better

        #And all the other objects and variables are updated:
        anc = anc_new
        phi_matrix = phi_matrix_new
        X = X_new
        p = p_new
        V = V_new
        inv_V = inv_V_new
        conditional = conditional_new
      }

      # Update beta whether tree accepted or not
      curr_trees[[j]] = simulate_beta_tilde(curr_trees[[j]],
                                      X,
                                      current_partial_residuals,
                                      sigma2,
                                      inv_V,
                                      tau_b,
                                      nu,
                                      ancestors)

      # Obtain the estimated beta's and subsequently the current tree fit
      beta_hat[[j]] = get_beta_hat(curr_trees[[j]])
      tau_vec[[j]] = tau[[j]]
      current_fit = X%*%beta_hat[[j]]
      pred_mat[,j] = current_fit

      # Calculate the new predictions and store the current fit
      predictions = predictions - tree_fits_store[,j] # subtract the old fit
      predictions = predictions + current_fit # add the new fit
      tree_fits_store[,j] = current_fit # update the new fit

    } # End loop through trees

    # Update sigma2 (variance of the residuals)
    sum_of_squares = sum((y_scale - predictions)^2)
    sigma2 = update_sigma2(sum_of_squares, n = length(y_scale), nu, lambda)

    # Update the variance of the terminal node parameters
    beta_trees = paste_betas(curr_trees,ntrees)
    tau_b = simulate_tau_b(beta_trees, sigma2, a,b)

    # Get the overall log likelihood
    log_lik = sum(dnorm(y_scale, mean = predictions, sd = sqrt(sigma2), log = TRUE))


  }# End iterations loop

  cat('\n') # Make sure progress bar ends on a new line

  return(list(trees = tree_store,
              sigma2 = sigma2_store*y_sd^2,
              y_hat = y_hat_store*y_sd + y_mean,
              beta_trees = beta_store,
              tau_trees = tau_store,
              tau_b_trees = tau_b_store,
              log_lik = log_lik_store,
              center_x = center,
              scale_x = scale,
              npost = npost,
              nburn = nburn,
              nthin = nthin,
              ntrees = ntrees,
              y_mean = y_mean,
              y_sd = y_sd,
              ancestors = ancestors
              ))
}

