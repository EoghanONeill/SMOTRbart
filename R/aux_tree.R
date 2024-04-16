# -------------------------------------------------------------------------#
# Description: this script contains auxiliar functions needed to update    #
# the trees with details and to map the predicted values to each obs       #
# -------------------------------------------------------------------------#

# 1. fill_tree_details: takes a tree matrix and returns the number of obs in each node in it and the indices of each observation in each terminal node
# 2. get_predictions: gets the predicted values from a current set of trees
# 3. get_children: it's a function that takes a node and, if the node is terminal, returns the node. If not, returns the children and calls the function again on the children
# 4. resample: an auxiliar function
# 5. get_ancestors: get the ancestors of all terminal nodes in a tree
# 6. update_s: full conditional of the vector of splitting probability.
# 7. update_vars_intercepts_slopes: updates the variances of the intercepts and slopes

# Fill_tree_details -------------------------------------------------------

fill_tree_details = function(curr_tree, X) {

  # Collect right bits of tree
  tree_matrix = curr_tree$tree_matrix

  # Create a new tree matrix to overwrite
  new_tree_matrix = tree_matrix

  # Start with dummy node indices
  node_indices = rep(1, nrow(X))

  # For all but the top row, find the number of observations falling into each one
  for(i in 2:nrow(tree_matrix)) {

    # Get the parent
    curr_parent = as.numeric(tree_matrix[i,'parent'])

    # Find the split variable and value of the parent
    split_var = as.numeric(tree_matrix[curr_parent,'split_variable'])
    split_val = as.numeric(tree_matrix[curr_parent, 'split_value'])

    # Find whether it's a left or right terminal node
    left_or_right = ifelse(tree_matrix[curr_parent,'child_left'] == i,
                           'left', 'right')
    if(left_or_right == 'left') {
      # If left use less than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] < split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] < split_val] = i
    } else {
      # If right use greater than condition
      new_tree_matrix[i,'node_size'] = sum(X[node_indices == curr_parent,split_var] >= split_val)
      node_indices[node_indices == curr_parent][X[node_indices == curr_parent,split_var] >= split_val] = i
    }
  } # End of loop through table

  return(list(tree_matrix = new_tree_matrix,
              node_indices = node_indices))

} # End of function

# Get predictions ---------------------------------------------------------

get_predictions = function(trees, X, single_tree = FALSE, ancestors) {

  # Stop nesting problems in case of multiple trees
  if(is.null(names(trees)) & (length(trees) == 1)) trees = trees[[1]]

  # Normally trees will be a list of lists but just in case
  if(single_tree) {

    # Deal with just a single tree
    if(nrow(trees$tree_matrix) == 1) {
      # predictions = rep(trees$tree_matrix[1, 'mu'], nrow(X))
      beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[1, 'beta_hat'],",")))
      predictions = rep(beta_hat[1], nrow(X))

    } else {
      # Loop through the node indices to get predictions
      predictions = rep(NA, nrow(X))
      unique_node_indices = unique(trees$node_indices)
      # Get the node indices for the current X matrix
      curr_X_node_indices = fill_tree_details(trees, X)$node_indices
      which_internal = which(trees$tree_matrix[,'terminal'] == 0)
      split_vars_tree <- trees$tree_matrix[which_internal, 'split_variable']

      if (ancestors == FALSE) {lm_vars <- c(1, sort(unique(as.numeric(split_vars_tree))))}
      #if (ancestors == 'all covariates') {lm_vars <- 1:ncol(X)}
      if (ancestors == TRUE) {get_ancs <- get_ancestors(trees)}

      n = nrow(X)

      # Now loop through all node indices to fill in details
      for(i in 1:length(unique_node_indices)) {
        if (ancestors == TRUE) {
          lm_vars = c(1, get_ancs[which(get_ancs[,'terminal'] == unique_node_indices[i]), 'ancestor']) # Get the corresponding ancestors of the current terminal node
        }
        X_node = matrix(X[,lm_vars], nrow=n)[curr_X_node_indices == unique_node_indices[i],]
        beta_hat = as.numeric(unlist(strsplit(trees$tree_matrix[unique_node_indices[i], 'beta_hat'],",")))
        predictions[curr_X_node_indices == unique_node_indices[i]] = X_node%*%beta_hat
      }
    }
    # More here to deal with more complicated trees - i.e. multiple trees
  } else {

    # Do a recursive call to the function
    partial_trees = trees
    partial_trees[[1]] = NULL # Blank out that element of the list
    predictions = get_predictions(trees[[1]], X, single_tree = TRUE, ancestors)  +
      get_predictions(partial_trees, X,
                      single_tree = length(partial_trees) == 1, ancestors)
    #single_tree = !is.null(names(partial_trees)))
    # The above only sets single_tree to if the names of the object is not null (i.e. is a list of lists)
  }

  return(predictions)
}

# get_children ------------------------------------------------------------

get_children = function(tree_mat, parent) {
  # Create a holder for the children
  all_children = NULL
  if(as.numeric(tree_mat[parent,'terminal']) == 1) {
    # If the node is terminal return the list so far
    return(c(all_children, parent))
  } else {
    # If not get the current children
    curr_child_left = as.numeric(tree_mat[parent, 'child_left'])
    curr_child_right = as.numeric(tree_mat[parent, 'child_right'])
    # Return the children and also the children of the children recursively
    return(c(all_children,
             get_children(tree_mat,curr_child_left),
             get_children(tree_mat,curr_child_right)))
  }
}

# Sample function ----------------------------------------------------------

resample <- function(x, ...) x[sample.int(length(x), size=1), ...]

# Get ancestors of each terminal node --------------------------------------

get_ancestors = function(tree){

  save_ancestor = NULL
  which_terminal = which(tree$tree_matrix[,'terminal'] == 1)

  if(nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          ancestor = NULL)
  } else {
    for (k in 1:length(which_terminal)){
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent'])) # get the 1st parent
      get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the covariate associated to the row of the parent

      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  # parent   = get_parent,
                                  ancestor = get_split_var))
      while (get_parent > 1){
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # then, get the subsequent parent
        get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the covariate associated to the row of the new parent
        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    # parent   = get_parent,
                                    ancestor = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }

  return(save_ancestor)
}

# update_s = function(var_count, p, alpha_s){
#   s_ = rdirichlet(1, alpha_s/p + var_count)
#   return(s_)
# }

update_s <- function(var_count, p, alpha_s) {
  # s_ = rdirichlet(1, as.vector((alpha_s / p ) + var_count))

  # // Get shape vector
  # shape_up = alpha_s / p
  shape_up = as.vector((alpha_s / p ) + var_count)

  # // Sample unnormalized s on the log scale
  templogs = rep(NA, p)
  for(i in 1:p) {
    templogs[i] = SoftBart:::rlgam(shape = shape_up[i])
  }

  if(any(templogs== -Inf)){
    print("alpha_s = ")
    print(alpha_s)
    print("var_count = ")
    print(var_count)
    print("templogs = ")
    print(templogs)
    stop('templogs == -Inf')
  }

  # // Normalize s on the log scale, then exponentiate
  # templogs = templogs - log_sum_exp(hypers.logs);
  max_log = max(templogs)
  templogs2 = templogs - (max_log + log(sum(exp( templogs  -  max_log ))))


  s_ = exp(templogs2)

  # if(any(s_==0)){
  #   print("templogs2 = ")
  #   print(templogs2)
  #   print("templogs = ")
  #   print(templogs)
  #   print("alpha_s = ")
  #   print(alpha_s)
  #   print("var_count = ")
  #   print(var_count)
  #   print("s_ = ")
  #   print(s_)
  #   stop('s_ == 0')
  # }

  ret_list <- list()
  ret_list[[1]] <- s_
  ret_list[[2]] <- mean(templogs2)


  return(ret_list)
}





update_alpha <- function(s, alpha_scale, alpha_a, alpha_b, p, mean_log_s) {

  # create inputs for likelihood

  # log_s <- log(s)
  # mean_log_s <- mean(log_s)
  # p <- length(s)
  # alpha_scale   # denoted by lambda_a in JRSSB paper

  rho_grid <- (1:1000)/1001

  alpha_grid <- alpha_scale * rho_grid / (1 - rho_grid )

  logliks <- alpha_grid * mean_log_s +
    lgamma(alpha_grid) -
    p*lgamma(alpha_grid/p) +
    (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)


  # logliks <- log(ddirichlet( t(matrix(s, p, 1000))  , t(matrix( rep(alpha_grid/p,p) , p , 1000)  ) ) ) +
  #   (alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid)
  # # dbeta(x = rho_grid, shape1 = alpha_a, shape2 = alpha_b, ncp = 0, log = TRUE)

  # logliks <- rep(NA, 1000)
  # for(i in 1:1000){
  #   logliks[i] <- log(ddirichlet(s  , rep(alpha_grid[i]/p,p) ) ) +
  #     (alpha_a - 1)*log(rho_grid[i]) + (alpha_b-1)*log(1- rho_grid[i])
  # }

  max_ll <- max(logliks)
  logsumexps <- max_ll + log(sum(exp( logliks  -  max_ll )))

  # print("logsumexps = ")
  # print(logsumexps)

  logliks <- exp(logliks - logsumexps)

  if(any(is.na(logliks))){
    print("logliks = ")
    print(logliks)

    print("logsumexps = ")
    print(logsumexps)

    print("mean_log_s = ")
    print(mean_log_s)

    print("lgamma(alpha_grid) = ")
    print(lgamma(alpha_grid))

    print("p*lgamma(alpha_grid/p) = ")
    print(p*lgamma(alpha_grid/p))

    print("(alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid) = ")
    print((alpha_a - 1)*log(rho_grid) + (alpha_b-1)*log(1- rho_grid))

    print("max_ll = ")
    print(max_ll)

    # print("s = ")
    # print(s)


  }

  # print("logliks = ")
  # print(logliks)

  rho_ind <- sample.int(1000,size = 1, prob = logliks)


  return(alpha_grid[rho_ind])
}





update_vars_intercepts_slopes <- function(trees, n_tress, sigma2, a0 = 1, b0 = 1, a1 = 1, b1 = 1){

  n_terminal = 0
  n_vars_terminal = 0
  sum_of_squares_inter = 0
  sum_of_squares_slopes = 0

  for (i in 1:n_tress) {
    # Get current set of trees
    tree = trees[[i]]
    # get the terminal nodes
    terminal_nodes = as.numeric(which(tree$tree_matrix[,'terminal'] == 1))
    # get all coefficients of the linear predictors for each terminal node
    all_coef = strsplit(tree$tree_matrix[terminal_nodes, 'beta_hat'], ',')
    # get intercepts
    inter = as.numeric(unlist(lapply(all_coef, '[', 1)))
    # get slopes
    slopes = as.numeric(unlist(lapply(all_coef, '[', -1)))

    n_terminal = n_terminal + length(terminal_nodes)
    n_vars_terminal = n_vars_terminal + length(slopes)
    sum_of_squares_inter = sum_of_squares_inter + sum(inter^2)
    sum_of_squares_slopes = sum_of_squares_slopes + sum(slopes^2)
  }
  return(list(var_inter = rgamma(1, (n_terminal/2) + a0, sum_of_squares_inter/(2*sigma2) + b0),
              var_slopes = rgamma(1, (n_vars_terminal/2) + a1, sum_of_squares_slopes/(2*sigma2) + b1)))
}


sample_move = function(curr_tree, i, nburn, trans_prob){

  if (nrow(curr_tree$tree_matrix) == 1 || i < max(floor(0.1*nburn), 10)) {
    type = 'grow'
  } else {
    type = sample(c('grow', 'prune', 'change'),  1, prob = trans_prob)
  }
  return(type)
}

# This function returns the Metropolis-Hastings acceptance probability in line
# with Kapelner, A., and Bleich, J. (2013). "bartMachine: Machine learning with
# Bayesian additive regression trees."

get_MH_probability <- function(curr_tree, new_tree,l_old, l_new,
                               # att_weights_current, att_weights_new,
                               # curr_partial_resid_rescaled,
                               # new_partial_resid_rescaled,
                               type, trans_prob,
                               alpha, beta) {

  # Number of terminal nodes in current tree
  b_j <- sum(curr_tree$tree_matrix[, "terminal"])

  # Get the tree type probabilities
  # if(nrow(curr_tree$tree_matrix) == 1 ){
  #   prob_grow <-  1 #trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }else{
  #   prob_grow <- trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }

  prob_grow <- trans_prob[1]
  prob_prune <- trans_prob[2]

  # print("type = ")
  # print(type)
  # Type is either "grow", "prune" or "change"
  if (type == "grow") {
    # The chosen node that it used to grow is the parent of the newly added last two rows in the NEW tree
    node_to_split <- as.numeric(new_tree$tree_matrix[nrow(new_tree$tree_matrix), "parent"])
    # Get the number of unique values for each covariate given the chosen node to split on in the CURRENT tree
    # n_unique <- get_n_unique(curr_tree, X, node_to_split)
    # Get the chosen covariate that it used to split on, which is the split variable in the NEW tree
    # covariate_to_split <- new_tree$tree_matrix[node_to_split, "split_variable"]
    # Obtain the additional parameters needed for the tree and transition ratios
    # p_j <- sum(n_unique >= 2 * node_min_size) # Number of available covariates to split on
    # n_j <- n_unique[covariate_to_split] - (node_min_size > 0) * (2 * node_min_size - 1) # Number of available values to split on given the chosen covariate, adjusted to the minimum leaf size
    w_2_star <- length(get_gen2(new_tree)) # Number of second generation nodes in the NEW tree

    # print("node_to_split = ")
    # print(node_to_split)
    #
    # print("curr_tree = ")
    # print(curr_tree)
    #
    # print("new_tree = ")
    # print(new_tree)

    depth_jl <- get_depth(curr_tree, node_to_split) # The depth of the chosen node to split on (in either the current or new tree)

    # tree_ratio <- alpha * (1 - alpha / ((2 + depth_jl)^beta))^2 / (((1 + depth_jl)^beta - alpha) * p_j * n_j)
    # transition_ratio <- prob_prune / prob_grow * b_j * p_j * n_j / w_2_star
    # print("depth_jl = ")
    # print(depth_jl)
    #
    # print("alpha * (1 - alpha / ((2 + depth_jl)^beta))^2 = ")
    # print(alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    #
    # print("((1 + depth_jl)^beta - alpha)  = ")
    # print(((1 + depth_jl)^beta - alpha) )
    #
    # print("alpha = ")
    # print(alpha)
    #
    # print("beta = ")
    # print(beta)

    tree_ratio <- alpha * (1 - (alpha / ((2 + depth_jl)^beta)) )^2 / (((1 + depth_jl)^beta - alpha) )
    transition_ratio <- (prob_prune / prob_grow) * (b_j  / w_2_star) # maybe this should be (w_2_star + 1) ?
  } else if (type == "prune") {
    # Finding the node used to prune is a bit more difficult then finding the node to split on
    pruned_node <- get_pruned_node(curr_tree, new_tree)
    # Get the number of unique values for each covariate given the chosen node it pruned on in the NEW tree
    # n_unique <- get_n_unique(new_tree, X, pruned_node)
    # Get the covariate that was used in the pruned node, which is the split variable in the CURRENT tree
    # pruned_covariate <- curr_tree$tree_matrix[pruned_node, "split_variable"]
    # Obtain the additional parameters needed for the tree and transition ratios
    # p_j_star <- sum(n_unique >= 2 * node_min_size) # Number of available covariates to split on
    # n_j_star <- n_unique[pruned_covariate] - (node_min_size > 0) * (2 * node_min_size - 1) # Number of available values to split on given the chosen covariate, adjusted to the minimum leaf size
    w_2 <- length(get_gen2(curr_tree)) # Number of second generation nodes in the CURRENT tree

    # print("pruned_node = ")
    # print(pruned_node)
    #
    # print("new_tree = ")
    # print(new_tree)

    depth_jl <- get_depth(new_tree, pruned_node) # The depth of the chosen node to split on (in either the current or new tree)

    # tree_ratio <- ((1 + depth_jl)^beta - alpha) * p_j_star * n_j_star / (alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    # transition_ratio <- prob_grow / prob_prune * w_2 / ((b_j - 1) * p_j_star * n_j_star)

    tree_ratio <- ((1 + depth_jl)^beta - alpha)  / (alpha * (1 - alpha / ((2 + depth_jl)^beta))^2)
    transition_ratio <- (prob_grow / prob_prune) * (w_2 / ((b_j - 1) ))
  } else {
    # For the changing step, the tree and transition ratios cancel out into a factor of 1
    transition_ratio <- 1
    tree_ratio <- 1
  }
  # l_new <- get_logL(new_tree, new_partial_resid_rescaled, att_weights_new, mu_mu, sigma2_mu, sigma2)
  # l_old <- get_logL(curr_tree, curr_partial_resid_rescaled, att_weights_current, mu_mu, sigma2_mu, sigma2)
  r <- exp(l_new - l_old) * transition_ratio * tree_ratio

  # print("l_new = ")
  # print(l_new)
  #
  # print("l_old = ")
  # print(l_old)
  #
  # print("transition_ratio = ")
  # print(transition_ratio)
  #
  # print("tree_ratio = ")
  # print(tree_ratio)

  return(min(1, r))
}



# This function returns the Metropolis-Hastings acceptance probability in line
# with Kapelner, A., and Bleich, J. (2013). "bartMachine: Machine learning with
# Bayesian additive regression trees."
get_MH_probability2 <- function(curr_tree, new_tree,l_old, l_new,
                                type, trans_prob,
                                alpha, beta) {
  # Number of terminal nodes in current tree
  # b_j <- sum(curr_tree$tree_matrix[, "terminal"])

  # Get the tree type probabilities
  # if(nrow(curr_tree$tree_matrix) == 1 ){
  #   prob_grow <-  1 #trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }else{
  #   prob_grow <- trans_prob[1]
  #   prob_prune <- trans_prob[2]
  # }

  prob_grow <- trans_prob[1]
  prob_prune <- trans_prob[2]
  # l_new <- get_logL(new_tree, new_partial_resid_rescaled, att_weights_new, mu_mu, sigma2_mu, sigma2)
  # l_old <- get_logL(curr_tree, curr_partial_resid_rescaled, att_weights_current, mu_mu, sigma2_mu, sigma2)

  if(type == 'grow'){
    # a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) )*
    #   ratio_grow(curr_tree, new_tree) * (prob_prune / prob_grow)
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log(ratio_grow(curr_tree, new_tree ) ))*
      (prob_prune / prob_grow)
  } else if(type == 'prune'){
    a = exp(l_new - l_old + get_tree_prior(new_tree, alpha, beta)  - get_tree_prior(curr_tree, alpha, beta) +
              log( ratio_prune(curr_tree, new_tree) ))*
      (prob_grow / prob_prune)
  } else{
    a = exp(l_new - l_old)
  }

  if(is.na(a)){
    print("get_tree_prior(new_tree, alpha, beta) = ")
    print(get_tree_prior(new_tree, alpha, beta))

    print("get_tree_prior(curr_tree, alpha, beta) ) = ")
    print(get_tree_prior(curr_tree, alpha, beta) )


    print("ratio_grow(curr_tree, new_tree ) ) = ")
    print(ratio_grow(curr_tree, new_tree ) )

    print("ratio_prune(curr_tree, new_tree)  = ")
    print(ratio_prune(curr_tree, new_tree) )

    print("new_tree = ")
    print(new_tree)

    print("curr_tree = ")
    print(curr_tree)

    print("alpha = ")
    print(alpha)

    print("beta = ")
    print(beta)

    print("l_new = ")
    print(l_new)

    print("l_old = ")
    print(l_old)


  }

  # r <- exp(l_new - l_old) * transition_ratio * tree_ratio
  return(min(1, a))
}

