#' modified log likelihood for Poisson distribution
#'
#' only calculate the non-zero elements of x likelihood
#' @param x is a vector
#' @param delta is a constant Poisson mean
#' @return loglikelihood of non-zero elements of x
dpoi <- function(x, delta)
{
  sum(dpois(x[x != 0], delta,log=T))
}


#' Main function for the nested EM algorithm for SMIXnorm.
#'
#' Performs the nested EM algorithm to find the MLE of the simplified mixture model.
#' @param dat input raw read count matrix. dim(dat)=J genes * I samples.
#' @param max_iter maximum number of iterations for the nested EM algorithm default is 20, recommend range (10, 50).
#' @param tol convergency criteria, default is 1e-2, recommend range (1e-5,1).
#' @return A list conatins the MLE of all parameters.
#'
#' @export
SMIXnorm.mle <- function(dat, max_iter=20, tol=1e-2)
{
  #return nonzero exit code is error occured
  exit.code <- 1
  cat('Start running SMIXnorm... \n')
  #check if the input arguments are correct
  if (sum(dat%%1!=0)!=0)
  {
    stop("Input data should be integer-valued read count matrix")
  }
  if (max_iter>50|max_iter<10)
  {
    stop("max_iter is out of range. Use proper maximum number of iteration (10, 50).")
  }
  if (tol>1|tol<1e-5)
  {
    stop("tol is out of range. Use proper convergence criteria (1e-5, 1)")
  }
  #end of checking
  #deal with NA
  dat[is.na(dat)] <- 0
  row.names(dat) <- 1:dim(dat)[1]
  #log-transformed data
  log_dat <- as.data.frame(log(dat+1))
  not_zero <- dat != 0
  #I = #of zero's in each row
  I <- rowSums(not_zero)
  samp_no <- dim(dat)[2]
  cat('Running iteration 1 ... \n')
  #set initial values for nested EM algorithm
  mu_i <- mu_i_final <- colMeans(log_dat)
  sigma_i <- sigma_i_final <- apply(log_dat, 2, sd)
  #proportion of 0's from point mass of each gene
  zip_pi <- pi_j <- 1-I/samp_no
  #proportion of genes that are truely expressed
  phi <- nonzero_prop <- mean(pi_j < .5)
  #poisson means for background read of genes not expressed
  delta <- poi_mean <- mean(colMeans(dat[pi_j > .5,]))
  zip_pi[I==0] <- pi_j[I==0] <- 1
  zip_pi[I==dim(dat)[2]] <- pi_j[I==dim(dat)[2]] <- 0
  #mu_i_final, sigma_i_final, phi, delta and pi_j are parameters after update
  #run 1 iteration of E step
  #outer E step: compute E(D_j|obs)
  log_numtor <- as.numeric(colSums(log(dnorm(t(log_dat), mu_i, sigma_i))))
  log_dentor1 <- as.vector(samp_no-I) * log(as.vector(pi_j) + (1-as.vector(pi_j)) * exp(-delta)) + I * log(1-as.vector(pi_j)) + as.vector(apply(dat, 1, dpoi, delta))
  #correct for numerical issues
  id_num <- which(log_numtor <= log_dentor1 & log_dentor1 <= (-700))
  log_numtor[id_num] <- log_numtor[id_num]-log_dentor1[id_num]
  log_dentor1[id_num] <- 0
  id_den <- which(log_dentor1 <= log_numtor & log_numtor <= (-700))
  log_dentor1[id_den] <- log_dentor1[id_den]-log_numtor[id_den]
  log_numtor[id_den] <- 0
  numtor <- phi*exp(log_numtor)
  dentor1 <- (1-phi)*exp(log_dentor1)
  dentor = numtor+dentor1
  D <- numtor/dentor
  D[is.na(D)] <- 1
  D[I==dim(dat)[2]] <- 1
  D[I==0] <- 0
  #end of outer E step
  #inner E step: compute E(Z_ij|obs,D_j)
  Z <- matrix(rep(pi_j^(1-D)/(pi_j^(1-D)+(1-pi_j)^(1-D)*exp(-delta*(1-D))),dim(dat)[2]),dim(dat)[1],dim(dat)[2])
  Z[dat!=0] <- 0
  #end of inner E step
  iter <- 1
  tolerance <- 9999
  while (tolerance > tol)
  {
    iter <- iter + 1
    if (iter > max_iter)
    {
      warning(paste('Algorithm exceeded maximum iteration (', max_iter, ') and stopped.'))
      break
    }
    cat(paste('Running iteration',iter,  '...\n'))
    #outer M step
    #update phi
    phi <- mean(D)
    #update mu_i and sigma_i
    mu_i_final <- colSums(log_dat * D)/sum(D)
    sigma_i_final <- sqrt(colSums(t(t(log_dat)-mu_i_final)^2 * D)/sum(D))
    #end of outer M step
    #Inner EM to update (pi_j, delta) K=5 times
    K <- 1
    while(K<=5)
    {
      #inner M step
      pi_j <- rowMeans(Z)
      pi_j[I==0] <- pi_j[I==0] <- 1
      pi_j[I==dim(dat)[2]] <- pi_j[I==dim(dat)[2]] <- 0

      deltadentor <-  (1-Z)*(1-D)
      deltanumtor <- (1-Z)*dat*(1-D)
      delta <- sum(deltanumtor)/sum(deltadentor)
      #inner E step
      if (K!=5)
      {
        Z <- matrix(rep(pi_j^(1-D)/(pi_j^(1-D)+(1-pi_j)^(1-D)*exp(-delta*(1-D))),dim(dat)[2]),dim(dat)[1],dim(dat)[2])
        Z[dat!=0] <- 0
      }
      K <- K+1
    }
    #End of Inner EM
    #Outer E step
    log_numtor <- as.numeric(colSums(log(dnorm(t(log_dat), mu_i, sigma_i))))
    log_dentor1 <- as.vector(samp_no-I) * log(as.vector(pi_j) + (1-as.vector(pi_j)) * exp(-delta)) + I * log(1-as.vector(pi_j)) + as.vector(apply(dat, 1, dpoi, delta))
    id_num <- which(log_numtor <= log_dentor1 & log_dentor1 <= (-700))
    log_numtor[id_num] <- log_numtor[id_num]-log_dentor1[id_num]
    log_dentor1[id_num] <- 0
    id_den <- which(log_dentor1 <= log_numtor & log_numtor <= (-700))
    log_dentor1[id_den] <- log_dentor1[id_den]-log_numtor[id_den]
    log_numtor[id_den] <- 0
    numtor <- phi*exp(log_numtor)
    dentor1 <- (1-phi)*exp(log_dentor1)
    dentor <- numtor+dentor1
    D <- numtor/dentor
    D[is.na(D)] <- 1
    D[I==dim(dat)[2]] <- 1
    D[I==0] <- 0
    #end of outer E step
    tolerance <- max(abs(mu_i_final-mu_i),abs(sigma_i_final-sigma_i),abs(nonzero_prop-phi))
    #update all parameters
    mu_i <- mu_i_final
    sigma_i <- sigma_i_final
    nonzero_prop <- phi
    poi_mean <- delta
    zip_pi <- pi_j
  }
  if (iter <= max_iter)
  {
    cat(paste0('\n','Algorithm converges at ', iter,'th iteration.\n' ))
    #successivefully obtained MLE, set exit.code to zero
    exit.code <- 0
  }
  return(list(mu_i_final = mu_i_final,
              sigma_i_final = sigma_i_final,
              D = D,
              phi = phi,
              delta = delta,
              pi_j = pi_j,
              exit.code = exit.code))
}


#' Produce SMIXnorm normalized expression matrix
#'
#' Calls the SMIXnorm.mle function to obtain MLE of the simplified mixture model,
#' then produces the normalized expression matrix.
#' @param dat input raw read count matrix. dim of dat = J genes * I samples.
#' @param max_iter maximum number of iterations for the nested EM algorithm default is 20, recommend range (10, 50).
#' @param tol convergency criteria, default is 1e-2, recommend range (1e-5,1).
#' @param appr binary True of False, indicates if the approximate version of normalization should be used.
#' @return A list contains the normalized expression matrix (MIX_normalized_log), proportion of expressed genes (phi) and probabilities of being expressed for all genes (D).
#' @export
SMIXnorm <- function(dat,max_iter = 20, tol = 1e-2,appr=T)
{
  SMIX_res <- SMIXnorm.mle(dat,max_iter=max_iter, tol=tol)
  subt <- matrix(rep(SMIX_res$mu_i_final,each=dim(dat)[1]),dim(dat)[1],dim(dat)[2])
  if (!appr)
  {
    SMIX_normalized_log <- as.data.frame(t(SMIX_res$D*(log(dat + 1)-subt)))
  }
  if (appr)
  {
    SMIX_normalized_log <- as.data.frame(t((SMIX_res$D>0.5)*(log(dat + 1)-subt)))
  }
  names(SMIX_res$D) <- rownames(dat)
  if (SMIX_res$exit.code==0)
  {
    cat('Normalization finished successfully \n')
  }
  return(list(SMIX_normalized_log = t(SMIX_normalized_log),phi = SMIX_res$phi, D = SMIX_res$D))
}


