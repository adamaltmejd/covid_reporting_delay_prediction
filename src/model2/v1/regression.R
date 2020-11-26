
#' Building the covariate matrix for the prior of the time point
#' 
#' @param  new_cases - (T x T) new cases each time point 
#'
build_X_time <- function(cases_remanin){
  
  n <- dim(cases_remanin)[1]
  m <- dim(cases_remanin)[2]
  X_ones <- matrix(1, nrow = n, ncol = m)
  X_diag <- matrix(0, nrow = n, ncol = m)
  diag(X_diag) <- 1
  X_na <- matrix(0, nrow = n, ncol = m)
  X_t <-  matrix(0, nrow = n, ncol = m)
  for(i in 1:(n-1)){
    X_t[i,i:n] <- 0:(n-i) 
    for(j in (i+1):(m)){
      X_na[i,j] <- 1*is.na(cases_remanin[i, j-1])
    }
    
  }
  
  return(list(X_ones = X_ones,
              X_diag = X_diag,
              X_na   = X_na,
              X_t    = X_t))
}