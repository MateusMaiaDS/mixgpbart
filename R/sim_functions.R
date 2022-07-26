# Simulating functions
#' @export
two_d_sim_two_node_gp_sum_each_tree <- function(n, # Number of observations
                                                mu1 =  c(-10,0,10), # Mean of the first terminal node
                                                mu2  = c(5,20,-15) , # Mean of the second terminal nodde
                                                nu = 0.1, # Getting the \nu parameter
                                                phi = 0.1, # Getting the \phi parameter
                                                tau = 10,
                                                unif_sample = FALSE,
                                                seed = NULL # Setting the seed
){
  #
  # Setting the seed 
  set.seed(seed)
  
  # Defining the kernel function structure
  omega_function <- function(x, x_star = NULL, nu, phi) {
    
    # Calculating the square matrix
    if (is.null(x_star)) {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x))^2)  +
        diag(1e-8, nrow = nrow(x))
    } else {
      kernel_matrix <- (nu^-1) * exp(-1 / (2 * phi^2) * as.matrix(stats::dist(x, x_star) )^2)
    }
    
    # Getting the kernel matrix
    return(kernel_matrix)
  }
  
  
  # Generating the x axis
  if(unif_sample){
    x <- expand.grid(lat = stats::runif(n = round(sqrt(n)),min = -10,max = 10),
                     lon = stats::runif(n = round(sqrt(n)),min = -10,max = 10))
  } else {
    x <- expand.grid(lat = seq(-10,10,length.out = round(sqrt(n))),
                     lon = seq(-10,10,length.out = round(sqrt(n))))
    
  }
  # Creating the true response
  y <- numeric(nrow(x))
  
  # Getting observation from different tree rules
  tree_one_n1 <- length(which(x[,1]< x[,2]))
  tree_one_n2 <- length(which(x[,1]>= x[,2]))
  
  tree_two_n1 <- length(which(x[,1] < -x[,2]))
  tree_two_n2 <- length(which(x[,1] >= - x[,2]))
  
  tree_three_n1 <- length(which(x[,1] < 0))
  tree_three_n2 <- length(which(x[,1] >= 0))
  
  
  # ==== Sampling for the first treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_one_node_one <- omega_function(x = x[x[,1] < x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) 
  # Defining the mean
  mean_tree_one_node_one <- rep(mu1[1],tree_one_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < x[,2] ] <-  y [ x[,1] < x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_one,
                                                              Sigma = var_tree_one_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_one_node_two <-  omega_function(x = x[x[,1] >= x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi)
  # Defining the mean
  mean_tree_one_node_two <- rep(mu2[1],tree_one_n2)
  
  y [ x[,1] >= x[,2] ] <-  y [ x[,1] >= x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_one_node_two,
                                                                Sigma = var_tree_one_node_two)
  
  
  # ==== Sampling for the second treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_two_node_one <-  omega_function(x = x[x[,1] < - x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi)
  # Defining the mean
  mean_tree_two_node_one <- rep(mu1[2],tree_two_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < -x[,2] ] <-  y [ x[,1] < -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_one,
                                                                Sigma = var_tree_two_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_two_node_two <-  omega_function(x = x[x[,1] >= -x[,2], ,drop =FALSE],
                                                          nu = nu, 
                                                          phi = phi) 
  mean_tree_two_node_two <- rep(mu2[2],tree_two_n2)
  
  y [ x[,1] >= -x[,2] ] <-  y [ x[,1] >= -x[,2] ] + MASS::mvrnorm(n = 1,mu = mean_tree_two_node_two,
                                                                  Sigma = var_tree_two_node_two)
  
  
  # ==== Sampling for the third treee ======
  
  # -- Sampling for the first node 
  
  # Defining the variance
  var_tree_three_node_one <-  omega_function(x = x[x[,1] < 0, ,drop =FALSE],
                                                            nu = nu, 
                                                            phi = phi) 
  # Defining the mean
  mean_tree_three_node_one <- rep(mu1[3],tree_three_n1)
  
  
  # Sampling from the first node
  y [ x[,1] < 0 ] <-  y [ x[,1] < 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_one,
                                                      Sigma = var_tree_three_node_one)
  
  
  
  # -- Sampling for the second node
  
  # Defining the variance
  var_tree_three_node_two <-  omega_function(x = x[x[,1] >= 0, ,drop =FALSE],
                                                            nu = nu, 
                                                            phi = phi) 
  # Defining the mean
  mean_tree_three_node_two <- rep(mu2[3],tree_three_n2)
  
  y [ x[,1] >= 0 ] <-  y [ x[,1] >= 0 ] + MASS::mvrnorm(n = 1,mu = mean_tree_three_node_two,
                                                        Sigma = var_tree_three_node_two)
  
  y_true <-y 
  y <- y + stats::rnorm(n = length(y),mean = 0,sd = tau^(-1/2))
  
  return(as.matrix(data.frame(x, y,y_true)))
  
}

# New simulation scenario
sim_2d_two_node <- function(n, seed, sd_value){

  
  # Setting a seed
  
  # Creating x and y
  x <- replicate(stats::runif(n = n,min = -1,max = 1),n = 2)
  colnames(x) <- c("lon","lat")
  y <- numeric(n)
  
  # Getting the left and right nodes
  left_index <- which(x[,"lon"]<0)
  right_index <- which(x[,"lon"]>=0)
  
  y[left_index] <- exp( (-0.5*(x[left_index,"lon"]^2+x[left_index,"lat"]^2)))
  y[right_index] <- -exp( (-0.5*(x[right_index,"lon"]^2+x[right_index,"lat"]^2)))
  f <- y
  y <- f + stats::rnorm(sd = sd_value,n = n)
  
  
  sim_data <- cbind(x,y,f)
  return(sim_data)  
}

# spatial_model <- sim_2d_two_node(n = 100,seed = 42,sd_value = 0.1)
# x <- spatial_model[,1:2,drop = FALSE]
# y <- spatial_model[,3,drop = FALSE]
# 
# # Test for spatial model
# spatial_model_test <- sim_2d_two_node(n = 200,seed = 43,sd_value = 0.1)
# x_test <- spatial_model_test[,1:2,drop = FALSE]
# y_test <- spatial_model_test[,3,drop = FALSE]
# 
# 
# bart_mod <- dbarts::bart(x.train = x,y.train = y,keeptrees = TRUE)
# bart_pred <- predict(bart_mod,newdata = x_test)
# 
# bart_mod$sigma %>% hist
# 
# soft_bart <- SoftBart::softbart(X = x,Y = y,X_test = x_test)
# soft_bart$sigma %>% hist()
# 
# gpbart_mod <- gp_bart(x = x,y = y,number_trees = 10)
# gpbart_pred <- predict(gpbart_mod,x_test = x_test)
# 
# gpbart_pred$out$sd %>% hist
# 
# crps(y = y,means = colMeans(gpbart_pred$out$pred),sds = rep(mean(gpbart_pred$out$sd),100))$CRPS
# crps(y = y,means = colMeans(bart_mod$yhat.train),sds = rep(mean(bart_mod$sigma),100))$CRPS
# crps(y = y,means = colMeans(soft_bart$y_hat_train),sds = rep(mean(soft_bart$sigma),100))$CRPS
# 
# rmse(obs = y_test,pred = colMeans(bart_pred))
# rmse(obs = y_test,pred = colMeans(gpbart_pred$out$pred))
# rmse(obs = y_test,pred = colMeans(soft_bart$y_hat_test))
# 
# plot3d(cbind(x_test,y_test))


#=========== Reading all arguments to debugg ===============#
# number_trees = 2 # Setting the number of trees
# node_min_size = 15 # Min node size
# mu = 0
# alpha = 0.5 # Alpha from prior
# beta = 5 # Beta from prior
# tau = 1 # Tau from prior
# n_iter = 2000 # Number of iterations
# burn = 500 # Number of burn
# thin = 1 # Number of thin
# rotation = TRUE # If rotated lon and lat will be used in tree building
# theta = NULL # If theta is NULL then the rotation angle will be randomly selected
# seed = NULL # Alpha vector values from the Dirichlet prior
# scale_boolean = TRUE
# # This will be defining the nu the default value
# nu_vector = NULL
# a_tau = 3 # Prior from a_v_ratio gamma
# d_tau = 1 # Prior from d_v_ratio gamma
# discrete_phi_boolean = FALSE
# x_scale =  TRUE
# gp_variables = colnames(x_train)   # Selecting the GP-Variables
# K_bart = 2
# prob_tau = 0.9
# kappa = 0.5
# bart_boolean = TRUE 
# bart_number_iter = 250


