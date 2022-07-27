## GP-Bart
#' @useDynLib mixgpbart
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp RcppEigen
# ==================================#
# Objects to test the tree_complete_conditional function
# ==================================#

tree_complete_conditional_gpbart <- function(tree, residuals,
                                             nu , phi ,
                                             tau_mu,kappa) {

  # Selecting the terminal nodes
  terminal_nodes <- tree[names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- vapply(terminal_nodes, function(x) {
    length(x$train_observations_index)
  }, numeric(1))

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$train_observations_index]
  })

  log_likedensity <- mapply(residuals_terminal_nodes,
                            terminal_nodes, FUN = function(res,node){
    mvnfast::dmvn(X = res,
                          mu = matrix(0, nrow = length(res)),
                     sigma = (1-kappa)*node$Omega_matrix+
                             diag(1/tau,nrow = nrow(node$Omega_matrix))+
                             kappa/tau_mu)
  })


  return(list(log_posterior = sum(log_likedensity)))
}

# Updating the predictions for the residuals
update_residuals <- function(tree,
                             x_train,
                             x_test,tau_mu,kappa,
                             nu, phi, residuals, tau) {

  # New g (new vector prediction for g)
  residuals_train_new <- rep(NA, nrow(x_train))
  residuals_test_new <- rep(NA,nrow(x_test))

  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Residuals terminal nodes
  residuals_terminal_nodes <- lapply(terminal_nodes, function(x) {
    residuals[x$train_observations_index]
  })


  residuals_sample <- mapply(terminal_nodes,
                                   residuals_terminal_nodes,
                                   FUN = function(node,resid_val)
                                   {gp_main_mix(x_train = x_train[node$train_observations_index,,drop = FALSE],
                                                x_test = x_test[node$test_observations_index,,drop = FALSE],
                                                omega = node$Omega_matrix,
                                                res_vec = resid_val,
                                                phi = phi,
                                                nu = nu,
                                                tau = tau,
                                                tau_mu = tau_mu,
                                                kappa = kappa)},SIMPLIFY = FALSE)

  # Adding the mu values calculated
  for(i in seq_along(terminal_nodes)) {
    # Saving g
    residuals_train_new[terminal_nodes[[i]]$train_observations_index] <- residuals_sample[[i]]$train_sample
    residuals_test_new[terminal_nodes[[i]]$test_observations_index] <- residuals_sample[[i]]$test_sample

  }
  return(list(residuals_train = residuals_train_new,
              residuals_test = residuals_test_new))
}


# ==============#
# rBart-GP FUNCTION
#' Fit a GP-BART model
#'
#' @param x_train A named train matrix x of the covariates from the training data.
#' @param y A named matrix y of the response from the training data
#' @param x_test A named test matrix
#' @param number_trees Number of trees used in the GP-BART model.
#' @param node_min_size Minimum number of observations in a terminal node
#' @param mu Initial values for the \eqn{\mu} parameter
#' @param alpha Parameter \eqn{\alpha} from tree prior
#' @param beta Parameter \eqn{\beta} from tree prior
#' @param tau Initial value for \eqn{\tau} parameter
#' @param n_iter Total number of iterations in GP-BART
#' @param burn Number of burn MCMC iterations
#' @param thin Number of thinning from MCMC
#' @param rotation Boolean to decide if rotated splits will be used
#' @param theta Fixed value for \eqn{\theta}, if NULL a random value for \eqn{\theta} will be used.
#' @param seed Setting a seed for the model
#' @param scale_boolean Boolean to scale or not the response \eqn{y}
#' @param a_tau Scale parameter from \eqn{\tau} prior
#' @param d_tau Rate parameter from \eqn{\tau} prior
#' @param discrete_phi_boolean Boolean to decide if it will be used a discrete grid for \eqn{\phi} proposals
#' @param x_scale Boolean to scale x or not
#' @param nu_vector A constant vector of length of number of trees with the \eqn{\nu} constant value
#' @param gp_variables Covariates used to build the covariance matrix from GP
#' @param K_bart Prior parameter from \eqn{\tau_{\mu}}
#' @param prob_tau Prior parameter from \eqn{\tau}
#' @param kappa Weight parameter \eqn{\kappa} from the proposed mixture of BART and GP-BART
#' @param bart_boolean Boolean to decide or not to warmup the model with BART samples
#' @param bart_number_iter Number of BART initial iterations to be used
#'
#' @return A 'gpbart_GPBART' model object,
#' @export
#'
gp_bart <- function(x_train, y, x_test,
                    number_trees = 2, # Setting the number of trees
                    node_min_size = 15, # Min node size,
                    mu = 0,
                    alpha = 0.5, # Alpha from prior
                    beta = 5, # Beta from prior
                    tau = 1, # Tau from prior,
                    n_iter = 2000, # Number of iterations
                    burn = 500, # Number of burn
                    thin = 1, # Number of thin
                    rotation = TRUE, # If rotated lon and lat will be used in tree building
                    theta = NULL, # If theta is NULL, then the rotation angle will be randomly selected
                    seed = NULL, # Alpha vector values from the Dirichlet prior
                    scale_boolean = TRUE,
                    # This will be defining the nu the default value
                    nu_vector = NULL,
                    a_tau = 3, # Prior from a_v_ratio gamma
                    d_tau = 1, # Prior from d_v_ratio gamma,
                    discrete_phi_boolean = FALSE,
                    x_scale =  TRUE,
                    gp_variables = colnames(x_train),   # Selecting the GP-Variables
                    K_bart = 2,
                    prob_tau = 0.9,
                    kappa = 0.5, bart_boolean = TRUE, bart_number_iter = 250) {

  # Changing the node_min_size
  if(node_min_size>=nrow(x_train)){
    stop("Node min size is greater than the number of observation")
  }

  # This parameter is a "scale paramter" to the GP
  phi_vector = rep(0.1 / (sqrt(number_trees)), number_trees)


  # Creating the prediction elements to be stored
  predictions = matrix(0, nrow = number_trees, ncol = nrow(x_train))
  predictions_test = matrix(0, nrow = number_trees, ncol = nrow(x_test))
  predictions_list =  NULL

  # If there's only one covariate
  rotation <- ncol(x_train) != 1

  mean_x <- NULL
  sd_x <- NULL

  # Scaling the x
  if(x_scale){
    x_train_original <- x_train
    x_test_original <- x_test

    # Scaled version
    xscale <- scale(x_train)

    mean_x <- attr(xscale,"scaled:center")
    sd_x <- attr(xscale,"scaled:scale")

    xscale_train <- scale(x_train,center = mean_x,scale = sd_x)
    xscale_test <- scale(x_test, center = mean_x, scale = sd_x)

    x_train <- as.matrix(x_train)
    x_test <- as.matrix(x_test)

  } else{
    x_train_original <- x_train
    x_test_original <- x_test
  }

  # Adjusting the kappa (avoiding the Infinity error)
  if(kappa == 1 ){
    kappa <- kappa - 2*.Machine$double.eps
  }

  if(kappa == 0 ){
    kappa <- kappa + 2*.Machine$double.eps
  }

  # Getting the maximum and minimum values from a distance matrix
  distance_matrix_x <- symm_distance_matrix(m1 = x_train[,gp_variables, drop = FALSE])
  distance_range <- range(distance_matrix_x[upper.tri(distance_matrix_x)])
  distance_min <- sqrt(distance_range[1])
  distance_max <- sqrt(distance_range[2])

  # Setting seed
  set.seed(seed)
  acc_ratio <- 0
  acc_ratio_phi <- 0
  acc_ratio_nu <- 0

  # Saving a_min and b_max
  a_min <- NULL
  b_max <- NULL

  a_min <- min(y)
  b_max <- max(y)

  # Scale values
  if(scale_boolean) {
    # Normalizing y
    y_scale <- normalize_bart(y = y)

    # Defing the nu vector if not in default
    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4 * number_trees * K_bart^2), number_trees)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees * K_bart^2)
    tau_mu_gpbart <- tau_mu_bart

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x_train,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)

  } else {

    # Not scaling the y
    y_scale <- y

    if(is.null(nu_vector)) {
      # Defining the nu value values on the maximum and minimum
      nu_vector <- rep((4 * number_trees * K_bart^2)/((max(y_scale) - min(y_scale))^2), number_trees)
    } else if(length(nu_vector) == 1) {
      nu_vector <- rep(nu_vector, number_trees)
    }

    # Calculating \tau_{\mu} based on the scale of y
    tau_mu_bart <- (4 * number_trees * K_bart^2)/((max(y_scale) - min(y_scale))^2)
    tau_mu_gpbart <- tau_mu_bart

    # Getting the optimal tau values
    d_tau <- rate_tau(x = x_train,
                      y = y_scale,
                      prob = prob_tau,
                      shape = a_tau)
  }


  # Creating the likelihood object list
  likelihood_object <- list()

  # Getting the number of observations
  n <- length(y)

  # Creating the stor of accepted tree verbs and which split it was
  verb_store_list <- list()

  # Getting the current trees
  current_trees <- list()
  current_trees_proposal <- list()

  # Fixed trees with true parameters
  fixed_trees <- list()

  # Creating the current_partial_residuals
  current_partial_residuals_matrix <-
    current_predictions_matrix <- matrix(NA, nrow = number_trees, ncol = nrow(x_train))

  # Creating the predictions saving
  current_partial_residuals_list <- list()
  current_predictions_list <- list()

  # Error of the matrix
  if(is.null(colnames(x_train))) {
    stop("Insert a valid NAMED matrix")
  }

  if(node_min_size == 0) {
    stop("Node Minimum Size need to be greater than 0")
  }

  if(length(nu_vector) != number_trees) {
    stop("Insert a valid \\nu vector for the number of trees")
  }

  if(length(phi_vector) != number_trees) {
    stop("Insert a valid \\phi vector for the number of trees")
  }

  # Recommendation about the min_node_size
  if(node_min_size < 15) {
    warning("\n It is recommended that the node_min_size should be of at least 15 observations.", immediate.=TRUE)
  }

  # Storage containers
  store_size <- (n_iter - burn)
  tree_store <- vector("list", store_size)
  tau_store <- c()

  y_hat_store <-
    y_hat_store_proposal <- matrix(NA, ncol = length(y), nrow = store_size)


  y_hat_test_store <- matrix(NA, ncol = nrow(x_test), nrow = store_size)

  # Storing the likelihoods
  log_lik_store <-
    log_lik_store_fixed_tree <- rep(NA,store_size)

  # Getting the numbers
  loglike_fixed_tree_residuals <- numeric()
  loglike_tree_residuals <- numeric()
  loglike_fixed_tree_residuals_matrix <-
    loglike_tree_residuals_matrix <- matrix(NA, nrow = store_size, ncol = number_trees)

  full_cond_store <-
    phi_store <-
    phi_proposal_store <- matrix(NA, ncol = number_trees, nrow = store_size)
  phi_vector_proposal <- rep(0.1, number_trees)

  # Creating the list of trees stumps
  for(i in seq_len(number_trees)) {
    # Creating the fixed two split trees
    current_trees[[i]] <- stump(x_train = x_train,
                                x_test= x_test, mu = mu)
    current_trees_proposal[[i]] <- stump(x_train = x_train,
                                         x_test= x_test, mu = mu)
  }

  names(current_trees) <-
    names(current_trees_proposal) <- vapply(seq_len(number_trees), function(x) paste0("tree_", x), character(1)) # Naming each tree

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_iter,
    style = 3, width = 50,
    label = "Running GP-Sum-Sampler..."
  )


  # Setting initial values for phi vector
  for(i in seq_len(n_iter)) {

    utils::setTxtProgressBar(progress_bar, i)

    # Changing the bart boolean, when reach the maximum
    if(i >= bart_number_iter){
      bart_boolean <- FALSE
      tau_mu <- tau_mu_gpbart
    } else tau_mu <- tau_mu_bart

    if((i > burn) && ((i %% thin) == 0)) {

      # Saving the store of the other ones
      curr <- (i - burn) / thin
      tree_store[[curr]] <- lapply(current_trees, remove_omega_plus_I_inv)

      tau_store[curr] <- if(scale_boolean){
        tau/((b_max-a_min)^2)
      } else {
        tau
      }

      # Getting the posterior for y_hat_train
      y_hat_store[curr, ] <- if(scale_boolean){

        # Getting the unnormalized version from tau
        unnormalize_bart(colSums(predictions),a = a_min,b = b_max)
      } else {

        colSums(predictions)
      }


      # Getting the posterior for y_hat_test
      y_hat_test_store[curr,] <- if(scale_boolean){
        unnormalize_bart(colSums(predictions_test), a = a_min, b = b_max)
      } else {
        colSums(predictions_test)
      }

      # Saving the current partial
      current_partial_residuals_list[[curr]] <- current_partial_residuals_matrix

      # Saving the predictions
      current_predictions_list[[curr]] <- if(scale_boolean){
        unnormalize_bart(predictions, a = a_min, b = b_max)
      } else {
        predictions
      }

      phi_store[curr, ] <- phi_vector
      phi_proposal_store[curr, ] <- phi_vector_proposal
      verb_store_list[[curr]] <- verb_store
    }

    # Creating a boolean to create the first trees only using BART model
    if(bart_boolean) {

      # Verb Store
      verb_store <- data.frame(verb = rep(NA, number_trees),
                               accepted = rep(NA, number_trees),
                               identical = rep(NA, number_trees))

      for(j in seq_len(number_trees)) {

        # Getting the verb list

        # Calculating the residuals for each tree
        if(number_trees > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change", "swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- sample(c("grow", "grow_projection"),
                                                                                               size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate

        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x_train = x_train,
          x_test = x_test,
          gp_variables = gp_variables,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )

        new_trees[[j]] <- current_trees[[j]]

        # Checking if the update tree generated a valid tree, if not skip likelihood calculations
        if( !identical(current_trees[[j]],new_trees[[j]]) ){


          # Calculating the likelihood of the new tree
          likelihood_new <- tree_complete_conditional_bart(
            tree = new_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x_train = x_train, tau_mu = tau_mu, tau = tau,
          )

          # Calculating the likelihood of the old tree
          likelihood_old <- tree_complete_conditional_bart(
            tree = current_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x_train = x_train, tau_mu = tau_mu, tau = tau
          )

          # Extracting only the likelihood
          l_new <- likelihood_new +
            tree_prior(
              tree = new_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )

          # Extracting only the likelihood
          l_old <- likelihood_old +
            tree_prior(
              tree = current_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )

          # Getting the log of transitin prob
          log_transition <- log_transition_prob(current_tree = current_trees[[j]],
                                                new_tree = new_trees[[j]],verb = verb)

          # (log) Probability of accept the new proposed tree
          acceptance <- (l_new - l_old + log_transition)

          # If Storage or not based on thin and burn parameters
          if((i > burn) && ((i %% thin) == 0)) {
            full_cond_store[curr, j] <- l_old
          }

        } else { # Accepting exact same trees
          # Calculating the likelihood of the new and old tree
          likelihood_new <- likelihood_old <- tree_complete_conditional_bart(
            tree = new_trees[[j]], # Calculate the full conditional
            residuals_values = current_partial_residuals,
            x_train = x_train, tau_mu = tau_mu, tau = tau
          )
          acceptance <- 0
        }
        # In case of acceptance

        if(acceptance > 0 || acceptance > -stats::rexp(1)) {

          # Counting acc ratio
          acc_ratio <- acc_ratio + 1

          # Make changes if accept
          current_trees <- new_trees

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- TRUE


          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_new

        } else {

          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_old

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- FALSE

        } # End of accept if statement

        # To update the mu values
        current_trees[[j]] <- update_mu_bart(
          tree = current_trees[[j]],
          x_train = x_train,
          residuals = current_partial_residuals,
          tau = tau,
          tau_mu = tau_mu)

        # Updating the BART predictions
        bart_pred_aux <- update_predictions_bart(
          tree = current_trees[[j]], x_train = x_train,x_test = x_test
        )

        predictions[j, ] <- bart_pred_aux$pred_train
        predictions_test[j,] <- bart_pred_aux$pred_test

        # Updating the current residuals
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]

      } # End of Loop through the trees

      # =================
      # ATTENTION HERE!!!
      # =================
    } else { # Going over the case where the BART-boolean is no more valid

      # Verb Store
      verb_store <- data.frame(verb = rep(NA,number_trees),
                               accepted = rep(NA,number_trees),
                               identical = rep(NA,number_trees))

      for(j in seq_len(number_trees)) {

        # Calculating the residuals for each tree
        if(number_trees > 1) {

          # Calculating what Chipman called as R(j) = y - g_others_trees
          if(number_trees > 2) {
            current_partial_residuals <- y_scale - colSums(predictions[-j, , drop = FALSE])
          } else {
            current_partial_residuals <- y_scale - predictions[-j, ]
          }
        } else {
          current_partial_residuals <- y_scale
        }

        # Propose a new tree based on the verbs: grow/prune/change/swap
        if(rotation){
          verb <- sample(c("grow", "grow_projection", "prune", "change", "change_projection", "swap"),
                         prob = c(0.125,0.125,0.25,0.20,0.20,0.1), size = 1)
        } else{
          verb <- sample(c("grow", "prune", "change","swap"),
                         prob = c(0.25,0.25,0.4,0.1), size = 1)
        }

        # Case of rotation
        if(rotation){
          if(i < max(floor(0.1 * burn), 10) | length(current_trees[[j]]) == 1) verb <- sample(c("grow","grow_projection"),
                                                                                              size = 1) # Grow the tree for the first few iterations
        } else {
          if(i < max(floor(0.1 * burn), 10) || length(current_trees[[j]]) == 1) verb <- "grow"  # Grow the tree for the first few iterations
        }

        # GETTING A NEW TREE
        new_trees <- current_trees # Creating new trees to updated as candidate

        new_trees[[j]] <- update_tree_verb(
          tree = current_trees[[j]],
          x_train = x_train,
          x_test = x_test,
          gp_variables = gp_variables,
          node_min_size = node_min_size,
          verb = verb, rotation = rotation, theta = theta
        )


        # ==================== #
        # Getting the Omega Inverse the current and the future tree
        # ==================== #

        # Checking if the previous tree were identitcal or not
        if(!identical(current_trees[[j]],new_trees[[j]])){

          # Getting the inverse for the current terminal nodes
          current_trees[[j]] <- inverse_omega_plus_I(tree = current_trees[[j]],
                                                     x_train = x_train, tau = tau,
                                                     nu = nu_vector[j],
                                                     phi = phi_vector[j],
                                                     kappa = kappa,
                                                     tau_mu = tau_mu,
                                                     gp_variables = gp_variables,
                                                     number_trees = number_trees)

          # Getting the inverse for the new tree terminal nodes
          new_trees[[j]] <- inverse_omega_plus_I(tree = new_trees[[j]],
                                                 x_train = x_train,tau = tau,
                                                 nu = nu_vector[j],
                                                 phi = phi_vector[j],
                                                 kappa = kappa,tau_mu = tau_mu,
                                                 number_trees = number_trees,
                                                 gp_variables = gp_variables)

          # Calculating the likelihood of the new tree
          likelihood_new <- tree_complete_conditional_gpbart(
            tree = new_trees[[j]],  # Calculate the full conditional
            residuals = current_partial_residuals,
            tau_mu = tau_mu,
            nu = nu_vector[j], phi = phi_vector[j],
            kappa = kappa
          )

          # Calculating the likelihood of the old tree
          likelihood_old <- tree_complete_conditional_gpbart(
            tree = current_trees[[j]], # Calculate the full conditional
            residuals = current_partial_residuals,
            tau_mu = tau_mu,
            nu = nu_vector[j], phi = phi_vector[j],
            kappa = kappa
          )

          # Extracting only the likelihood
          l_new <- likelihood_new$log_posterior +
            tree_prior(
              tree = new_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )

          # Extracting only the likelihood
          l_old <- likelihood_old$log_posterior +
            tree_prior(
              tree = current_trees[[j]], # Calculate the tree prior
              alpha = alpha,
              beta = beta
            )

          # Getting the log of transitin prob
          log_transition <- log_transition_prob(current_tree = current_trees[[j]],
                                                new_tree = new_trees[[j]],verb = verb)

          # (log) Probability of accept the new proposed tree
          acceptance <- (l_new - l_old + log_transition)

          # If Storage or not based on thin and burn parameters
          if((i > burn) && ((i %% thin) == 0)) {
            full_cond_store[curr, j] <- l_old
          }
          # In case of them being identical
        } else {

          acceptance <- 0.1
          # Checking if none Omega matrix element were calculated
          if(is.null(current_trees[[j]][unlist(lapply(current_trees[[j]],
                                                      function(node) node$terminal==1))][[1]]$Omega_matrix)){
            # Creating the current tree object
            # Getting the inverse for the current terminal nodes
            current_trees[[j]] <- new_trees[[j]] <- inverse_omega_plus_I(tree = current_trees[[j]],
                                                                         x_train = x_train, tau = tau,
                                                                         nu = nu_vector[j],
                                                                         phi = phi_vector[j],
                                                                         kappa = kappa,
                                                                         tau_mu = tau_mu,
                                                                         number_trees = number_trees)

            # Creating the likelihood object
            likelihood_new <- likelihood_old <-
              tree_complete_conditional_gpbart(
                tree = current_trees[[j]], # Calculate the full conditional
                residuals = current_partial_residuals,
                tau_mu = tau_mu,
                nu = nu_vector[j], phi = phi_vector[j],
                kappa = kappa
              )
          } else{ # Replacing the likelihood_object
            likelihood_new <- likelihood_old <- likelihood_object[[j]]
          } # end of the if for null likelihood and current tree

        }

        if(acceptance > 0 || acceptance > -stats::rexp(1)) { #
          acc_ratio <- acc_ratio + 1

          # Checking whether the trees are identical
          if(identical(current_trees[[j]], new_trees[[j]])){
            verb_store[j,"identical"] <- TRUE
          } else {
            verb_store[j,"identical"] <- FALSE
          }

          # Make changes if accept
          current_trees <- new_trees

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- TRUE

          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_new

        } else {
          # Storing likelihood matrix objects
          likelihood_object[[j]] <- likelihood_old

          # Create a data.frame with the verb and if it was accepted or not
          verb_store[j,"verb"] <- verb
          verb_store[j,"accepted"] <- FALSE
          verb_store[j,"identical"] <- FALSE

        } # End of accept if statement


        # EQUATION FROM SECTION 4
        # ==== Using the prediction from R_star_bar

        update_residuals_aux <- update_residuals(
          tree = current_trees[[j]],
          x_train = x_train,x_test = x_test,
          residuals = current_partial_residuals,
          phi = phi_vector[j], nu = nu_vector[j], tau = tau,
          tau_mu = tau_mu,kappa = kappa
        )

        predictions[j, ] <- update_residuals_aux$residuals_train
        predictions_test[j, ]<- update_residuals_aux$residuals_test

        # # To update phi
        mh_update_phi <- update_phi_marginal(current_tree_iter = current_trees[[j]],
                                             residuals = current_partial_residuals,
                                             x_train = x_train,nu = nu_vector[j],
                                             phi = phi_vector[j],
                                             gp_variables = gp_variables,
                                             likelihood_object = likelihood_object[[j]],
                                             discrete_phi = discrete_phi_boolean,
                                             tau = tau,kappa = kappa,
                                             tau_mu = tau_mu,
                                             distance_min = distance_min,
                                             distance_max = distance_max)

        # In case of accept the update over \phi update everything
        if(mh_update_phi$phi_boolean) {

          # Updating the tree and the \phi object from the tree
          current_trees[[j]] <- mh_update_phi$tree

          # Updating the likelihood objects
          likelihood_object[[j]] <- mh_update_phi$likelihood_object

          # Updating the phi value
          phi_vector[j] <- mh_update_phi$phi_proposal

        } # If doesn't accept, nothing changes.



        # current_partial_residuals_matrix<-
        current_partial_residuals_matrix[j, ] <- current_partial_residuals
        current_predictions_matrix[j, ] <- predictions[j, ]

      } # End of Loop through the trees
    }

    tau <- update_tau_linero(x_train = x_train,
                             y = y_scale,
                             y_hat = colSums(predictions),
                             curr_tau = tau)

  } # End of Loop through the n_inter
  cat("\n")

  # Returning X to its original scale
  if(x_scale) {
    x_train <- x_train_original
    x_test <- x_test_original
  }

  results <- list(trees = tree_store,
                  tau_store = tau_store,
                  y_hat = y_hat_store,
                  y_hat_test = y_hat_test_store,
                  log_lik = log_lik_store,
                  log_lik_fixed_tree = log_lik_store_fixed_tree,
                  loglike_fixed_tree_residuals_matrix = loglike_fixed_tree_residuals_matrix,
                  loglike_tree_residuals_matrix = loglike_tree_residuals_matrix,
                  full_cond = full_cond_store,
                  phi_store = phi_store,
                  phi_proposal_store = phi_proposal_store,
                  nu_vector = nu_vector,
                  y = y_scale,
                  X = x_train_original,
                  x_scale = x_scale,
                  mean_x = mean_x,
                  sd_x = sd_x,
                  scale_boolean = scale_boolean,
                  acc_ratio = acc_ratio,
                  acc_ratio_phi = acc_ratio_phi,
                  iter = n_iter,
                  burn = burn,
                  thin = thin,
                  store_size = store_size,
                  number_trees = number_trees,
                  node_min_size = node_min_size,
                  a_min = a_min,
                  b_max = b_max,
                  a_tau = a_tau,
                  d_tau = d_tau,
                  current_partial_residuals_list = current_partial_residuals_list,
                  beta = beta,
                  current_predictions_list = current_predictions_list,
                  tau_mu = tau_mu, kappa = kappa,
                  verb_store_list = verb_store_list)
  class(results) <- "gpbart_GPBART"
  return(results)
}

# #Do a MH for PHI
update_phi_marginal <- function(x_train, current_tree_iter,residuals,
                                seed = NULL,
                                tau,
                                kappa,
                                tau_mu,
                                phi, nu,
                                likelihood_object, p, gp_variables,
                                discrete_phi = TRUE,
                                distance_min,
                                distance_max) {

  # Increased the range of tree proposal
  if(discrete_phi){
    phi_proposal <- sample(c(0.1,0.5,1,5,10), size = 1)
  } else {
    phi_proposal <- stats::runif(1, min = distance_min, max = distance_max)
  }

  # Calculating the likelihood from the new step
  tree_from_phi_proposal <- inverse_omega_plus_I(tree = current_tree_iter,
                                                 x_train = x_train,nu = nu, tau = tau,
                                                 kappa = kappa,tau_mu = tau_mu,
                                                 phi = phi_proposal, gp_variables = gp_variables)

  likelihood_phi_proposal <- tree_complete_conditional_gpbart(tree = tree_from_phi_proposal,
                                                              residuals = residuals,
                                                              nu = nu, tau_mu = tau_mu,kappa = kappa,
                                                              phi = phi_proposal)

  # Old phi likelhood
  l_old_phi <- likelihood_object$log_posterior

  # Proposal likelihood
  l_proposal_phi <- likelihood_phi_proposal$log_posterior

  # (log) Probability of accept the new proposed tree
  acceptance_phi <- l_proposal_phi - l_old_phi

  # If storage for phi

  if(acceptance_phi > 0 || acceptance_phi > -stats::rexp(1)) { #

    # Nu boolean to see if was accepted or not
    phi_boolean <- TRUE
    return(list(phi_boolean = phi_boolean,
                likelihood_object = likelihood_phi_proposal,
                tree = tree_from_phi_proposal,
                phi_proposal = phi_proposal)) # Returning the proposal value for phi
  } else {
    # Case of not accepting
    phi_boolean <- FALSE
    return(list(phi_boolean = phi_boolean)) # Returning the old value for phi
  } #
}


# Function to return the depth trees

tree_depth_hist <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for (i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- max(vapply(tree, "[[", numeric(1), "depth_node"))
    }
  }
  return(tree_depth)
}

# Get the covariate splits
tree_var_hist <- function(gpbart_model) {
  tree_depth <- c()

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for(i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth <- c(tree_depth, sapply(tree, "[[", "node_var"))
    }
  }
  return(tree_depth)
}

# Function to count th enumber of terminal nodes in a tree
tree_count_terminals <- function(gpbart_model) {
  tree_depth <- matrix(NA, nrow = length(gpbart_model$trees), ncol = gpbart_model$number_trees)

  for(k in seq_along(gpbart_model$trees)) {
    tree_iter <- gpbart_model$trees[[k]]
    for(i in seq_len(gpbart_model$number_trees)) {
      tree <- tree_iter[[i]]
      tree_depth[k, i] <- sum(vapply(tree, "[[", numeric(1), "terminal"))
    }
  }
  return(tree_depth)
}

# Tau values
get_tau_values <- function(gpbart_model) {

  # n_iter
  n_iter <- length(gpbart_model$trees)

  tau_iters <- rep(list(NA), gpbart_model$number_trees)

  num_trees <- gpbart_model$number_trees

  for(j in seq_len(n_iter)) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for tau terminals
    tau_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(tau_terminals) <- paste0("node_", 0:49)

    for(k in seq_len(num_trees)) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)

      for(nodes in seq_along(all_nodes)) {
        if(tree[[all_nodes[nodes]]]$terminal == 1) {
          # tree[[nodes]]$tau <- sample(c(0.001,seq(0.005,10,by = 0.005),10,50,100), 1)
          tau_terminals[k, as.numeric(stringr::str_extract(all_nodes[nodes], pattern = "\\(?[0-9,.]+\\)?")) + 1] <- tree[[all_nodes[nodes]]]$tau
        }
      }
    }

    # Saving all
    tau_iters[[j]] <- tau_terminals
  }
  return(tau_iters)
}

# mu values
get_mu_values <- function(gpbart_model) {

  # n_iter
  n_iter <- length(gpbart_model$trees)

  mu_iters <- rep(list(NA), gpbart_model$number_trees)

  num_trees <- gpbart_model$number_trees

  for(j in seq_len(n_iter)) {
    num_trees <- length(gpbart_model$trees[[j]])
    # Creating dummy for mu terminals
    mu_terminals <- matrix(NA, nrow = num_trees, ncol = 50)
    colnames(mu_terminals) <- paste0("node_", 0:49)

    for(k in seq_len(num_trees)) {
      # Tree
      tree <- gpbart_model$trees[[j]][[k]]
      all_nodes <- names(tree)

      for(nodes in seq_len(all_nodes)) {
        if (tree[[all_nodes[nodes]]]$terminal == 1) {
          # tree[[nodes]]$mu <- sample(c(0.001,seq(0.005,10,by = 0.005),10,50,100), 1)
          mu_terminals[k, as.numeric(stringr::str_extract(all_nodes[nodes], pattern = "\\(?[0-9,.]+\\)?")) + 1] <- tree[[all_nodes[nodes]]]$mu
        }
      }
    }

    # Saving all
    mu_iters[[j]] <- mu_terminals
  }
  return(mu_iters)
}


# Getting Omega Inverse + Diag Inverse
inverse_omega_plus_I <- function(tree,
                                 x_train ,
                                 nu, phi,kappa,
                                 tau,tau_mu,
                                 number_trees = number_trees,
                                 gp_variables = colnames(x_train)  # Selecting which gp-variables to use

) {
  # Selecting terminal nodes names
  names_terminal_nodes <- names(which(vapply(tree, "[[", numeric(1), "terminal") == 1))

  # Selecting the terminal nodes
  terminal_nodes <- tree[names_terminal_nodes]

  # Number of nodes
  n_node <- length(terminal_nodes)

  # Picking each node size
  nodes_size <- sapply(terminal_nodes, function(x) {
    length(x$train_observations_index)
  })

  # Calculating Omega matrix INVERSE
  distance_matrices <- mapply(terminal_nodes, FUN = function(y) {
    symm_distance_matrix(matrix(x_train[y$train_observations_index, gp_variables], nrow = length(y$train_observations_index)))
  }, SIMPLIFY = FALSE)

  # Calculating Omega
  Omega_matrix <- mapply(distance_matrices, FUN = function(dist_m) {
    kernel_function(
      squared_distance_matrix = dist_m,
      nu = nu, phi = phi)
  }, SIMPLIFY = FALSE)

  # Checking if diagonal
  is_Omega_diag <- lapply(Omega_matrix, is_diag_matrix)

  # Adding the Omega_matrix_plus_I_Inv
  for(i in seq_along(names_terminal_nodes)) {
    tree[[names_terminal_nodes[i]]]$Omega_matrix <- Omega_matrix[[names_terminal_nodes[i]]]
    tree[[names_terminal_nodes[i]]]$is_Omega_diag <- is_Omega_diag[[names_terminal_nodes[i]]]
  }
  return(tree)
}

# # Removing the Omega_plus_I_inv object
remove_omega_plus_I_inv <- function(current_tree_iter) {

  # Selecting terminal nodes names
  # names_terminal_nodes <- names(which(vapply(current_tree_iter, "[[", numeric(1), "terminal") == 1))
  names_terminal_nodes <- names(current_tree_iter)

  for(i in names_terminal_nodes) {
      current_tree_iter[[i]]$Omega_matrix <-
      current_tree_iter[[i]]$is_Omega_diag <- NULL
  }
  return(current_tree_iter)
}




# Get train predictions
get_train_predictions <- function(gpbart_mod) {

  # Getting the quantile
  gpbart_sum_pred <- do.call(rbind, lapply(gpbart_mod$current_predictions_list, colSums))

  # Returning the matrix of final predictions
  return(gpbart_sum_pred)
}

# Calculating a PI coverage
#' @export
pi_coverage <- function(y, y_hat_post, sd_post,only_post = FALSE, prob = 0.5,n_mcmc_replications = 1000){

  # Getting the number of posterior samples and columns, respect.
  np <- nrow(y_hat_post)
  nobs <- ncol(y_hat_post)

  full_post_draw <- list()

  # Setting the progress bar
  progress_bar <- utils::txtProgressBar(
    min = 1, max = n_mcmc_replications,
    style = 3, width = 50 )

  # Only post matrix
  if(only_post){
    post_draw <- y_hat_post
  } else {
    for(i in 1:n_mcmc_replications){
      utils::setTxtProgressBar(progress_bar, i)

      full_post_draw[[i]] <-(y_hat_post + replicate(sd_post,n = nobs)*matrix(stats::rnorm(n = np*nobs),
                                                                             nrow = np))
    }
  }

  if(!only_post){
    post_draw<- do.call(rbind,full_post_draw)
  }

  # CI boundaries
  low_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = prob/2)})
  up_ci <- apply(post_draw,2,function(x){stats::quantile(x,probs = 1-prob/2)})

  pi_cov <- sum((y<=up_ci) & (y>=low_ci))/length(y)

  return(pi_cov)
}
