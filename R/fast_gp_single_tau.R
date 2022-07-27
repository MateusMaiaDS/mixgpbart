# GP-function main
gp_main_slow <- function(x_train, y_train, x_star, tau,
                           phi, nu, distance_matrix_train,
                           get_sample =  FALSE) {


  # Getting the distance matrix from x_train and x_star
  distance_matrix_K_star <- distance_matrix(m1 = x_train, m2 = x_star)
  distance_matrix_K_star_star <- symm_distance_matrix(m1 = x_star)

  # Calculating the K elements from the covariance structure
  n_train <- nrow(x_train)
  if(tau < 1e13){
    K_y <- kernel_function(squared_distance_matrix = distance_matrix_train,
                           nu = nu,
                           phi = phi) + diag(x = 1/(tau), nrow = n_train)
  } else {
    K_y <- PD_chol(kernel_function(squared_distance_matrix = distance_matrix_train,
                           nu = nu,
                           phi = phi)) + diag(x = 1/(tau), nrow = n_train)
  }
  K_diag <- is_diag_matrix(K_y)
  K_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star,
                            nu = nu, phi = phi)

  # Calculating \alpha
  if(K_diag) {
    L <- diag(K_y)
    alpha <- y_train/L
  } else {
    L <- chol(K_y)
    alpha <- backsolve(L, backsolve(L, y_train, transpose = TRUE, k = n_train), k = n_train)
  }
  mu_star <- crossprod(K_star, alpha)

  # this line is fucking up everything
  # mu_star <- matrix(mu_star,nrow = n_train)
  # print(mu_star[1:5])

  # Here the abs is because the smallest values that are coming from here are due to numerical approximations.
  if(isTRUE(get_sample)) {

    K_star_star <- kernel_function(squared_distance_matrix = distance_matrix_K_star_star,
                                   nu = nu, phi = phi)

    cov_star <- K_star_star - crossprod(K_star,solve(K_y,K_star))

    residuals_sample <- rMVN_var(mean = mu_star,Sigma = cov_star)

    results <- list(mu_pred = unlist(residuals_sample))

  } else {
    results <- list(mu_pred = mu_star)
  }

  # ===============#
  return(results)
}

gp_main_mix <- function(x_train,
                        x_test,
                        res_vec,
                        omega,
                        kappa,
                        tau,
                        nu,
                        phi,
                        tau_mu){

  # Getting n_train 
  n_train <- nrow(omega)
  
  # Calculating the mean for the the training
  K_y <- diag(1/tau,nrow = nrow(omega)) + (1-kappa)*omega + kappa/tau_mu
  K_star_train <- (1-kappa)*omega + kappa/tau_mu
  
  L <- chol(K_y)
  alpha <- backsolve(L, backsolve(L, res_vec, transpose = TRUE, k = n_train), k = n_train)
  
  mu_star_train <- crossprod(K_star_train,alpha)
  
  # Getting the covariance matrix
  v <- backsolve(L, K_star_train, transpose = TRUE, k = n_train)
  
  cov_star_train <- ((1-kappa)*omega + kappa/tau_mu) + crossprod(v)
  
  
  train_sample <-rMVN_var(mean = mu_star_train,Sigma = cov_star_train)
  
  # Doing the same now, but for test predictions
  K_star_test <- (1-kappa)*kernel_function(squared_distance_matrix = distance_matrix(m1 = x_train,m2 = x_test),
                                 nu = nu,phi = phi) + kappa/tau_mu
  
  K_star_star <- (1-kappa)*kernel_function(squared_distance_matrix = symm_distance_matrix(m1 = x_test),
                                 nu = nu, phi = phi) + kappa/tau_mu
  
  mu_star_test <- crossprod(K_star_test,alpha)
  
  v_test <- backsolve(L,K_star_test,transpose = TRUE,k = n_train)
  
  cov_star_test <- K_star_star-crossprod(v_test)
  
  test_sample <- rMVN_var(mean = mu_star_test,Sigma = cov_star_test)
  

  return(list(train_sample = train_sample,
              test_sample = test_sample))

}

simple_gp <- function(x_train,
                      x_test,
                      train_sample,
                      nu,
                      phi){

  # Getting the mean for new observations
  omega <-  kernel_function(symm_distance_matrix(m1 = x_train),
                            nu = nu,phi = phi)

  omega_new <- kernel_function(squared_distance_matrix = distance_matrix(m1 = x_train,m2 = x_test),
                               nu = nu,phi = phi)
  # omega_new <- omega_new + matrix(1/tau_mu,nrow = nrow(omega_new), ncol = ncol(omega_new))

  L <- PD_chol(omega)
  # alpha <- backsolve(L, backsolve(L, matrix(train_sample,nrow = nrow(x_train)), transpose = TRUE, k = nrow(x_train)), k = nrow(x_train))
  alpha <- backsolve(L, backsolve(L, train_sample, transpose = TRUE, k = nrow(x_train)), k = nrow(x_train))

  # print(omega_new)
  mu_star <- crossprod(omega_new,solve(omega+diag(1e-9,nrow = nrow(omega)),train_sample))

  return(mu_star)
}

# Function to create the the function K that will be used
# in a Gaussian process (Andrew's Version)
kernel_function <- function(squared_distance_matrix, nu, phi) {

  # Calculating the square matrix
  kernel_matrix <- (exp(-squared_distance_matrix / (2 * phi^2))) / nu

  # Case nu = 0
  if(nu == 0 || nu > 1e13){
    kernel_matrix <- matrix(0, nrow = dim(squared_distance_matrix)[1],
                            ncol = dim(squared_distance_matrix)[2])
  }
  # Getting the kernel matrix
  return(kernel_matrix)
}
