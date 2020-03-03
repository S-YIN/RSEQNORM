#' log likelihood for zero inflated Poisson
#'
#' Input is a matrix with column specific Poisson mean and row specific probability of extra zero
#' @param x is a matrix, x[j,i] follows ZIP(delta[i],pi_j[j])
#' @param delta is a vector for poisson means
#' @param pi_j is a vector for prob. of extra zero
#' @return the log-likelihood of each row of x
dzip <- function(x, delta,pi_j)
{
  # zip_llh = ZIP log-likelihood matrix for I*J data points
  zip_llh <- matrix(NA,dim(x)[1],dim(x)[2])
  for (j in 1: dim(x)[1])
  {
    indictr <- ifelse(x[j,]>0,1,0)
    # if data contain zeros, use ZIP pmf
    if (sum(indictr)!=0)
    {
      zip_llh[j,] <- indictr*(log(1-pi_j[j])+dpois(as.numeric(x[j,]),delta,T))+(1-indictr)*(log(pi_j[j]+(1-pi_j[j])*(exp(-delta))))
    }
    # if no zero, use Poisson pmf
    if (sum(indictr)==0)
    {
      zip_llh[j,] <- (log(pi_j[j]+(1-pi_j[j])*(exp(-delta))))
    }
  }
  den_zip <- rowSums(zip_llh)
  den_zip
}


#' Main function for the nested EM algorithm for MIXnorm.
#'
#' Performs the nested EM algorithm to find the MLE of the mixture model.
#' @param dat input raw read count matrix. dim(dat)=J genes * I samples.
#' @param max_iter maximum number of iterations for the nested EM algorithm default is 20, recommend range (10, 50).
#' @param tol convergency criteria, default is 1e-2, recommend range (1e-5,1).
#' @return A list contains the MLE of all parameters.
#'
#' @export
#' @import truncnorm
MIXnorm.mle <- function(dat, max_iter=20, tol=1e-2)
{
  #return nonzero exit code is error occured
  exit.code <- 1
  cat('Start running MIXnorm... \n')
  #check if the input arguments are correct
  if (sum(dat%%1!=0)!=0)
  {
    stop("Input data should be integer-valued read count matrix")
  }

  if (max_iter>50|max_iter<10)
  {
    stop("max_iter is out of range. Use proper maximum number of iteration (10, 50).")
  }

  if (tol>1|max_iter<1e-5)
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
  delta <- poi_mean <- colMeans(dat[pi_j > .5,])
  delta[delta==0] <- 0.001
  zip_pi[I==0] <- pi_j[I==0] <- 1
  zip_pi[I==dim(dat)[2]] <- pi_j[I==dim(dat)[2]] <- 0
  #mu_i_final, sigma_i_final, phi, delta and pi_j are parameters after update
  #run 1 iteration of E step
  #outer E step: compute E(D_j|obs)
  log_numtor <- apply(matrix(log(truncnorm::dtruncnorm(t(log_dat),a=0,b=Inf, mu_i, sigma_i)),dim(log_dat)[2],dim(log_dat)[1]), 2, sum)
  log_dentor1 <- dzip(x=dat,delta=delta,pi_j=pi_j)
  #correct for numerical issues
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
  #inner E step: compute E(Z_ij|obs,D_j)
  Z <- dat
  for (i in 1:dim(dat)[2])
  {
    Z[,i] <- pi_j^(1-D)/(pi_j^(1-D)+(1-pi_j)^(1-D)*exp(-delta[i]*(1-D)))
  }
  Z[dat!=0] <- 0
  #inner E step: compute E(T_i|obs,D)
  beta <- -mu_i/sigma_i
  T_i <- sum(D)*pnorm(beta)/(1-pnorm(beta))
  #sum of augment data
  mis_sum <- T_i*(mu_i-dnorm(-beta)*sigma_i/pnorm(beta))
  #sum of augment data squares
  mis_sum_sq <- T_i*(sigma_i^2+mu_i^2+sigma_i^2*(-beta)*dnorm(-beta)/pnorm(beta)-2*mu_i*sigma_i*dnorm(beta)/pnorm(beta))
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
    #end of outer M step
    #Inner EM to update (pi_j, delta_i) and (mu_i,sigma_i) K=5 times
    K <- 1
    while(K<=5)
    {
      #inner M step
      pi_j <- rowMeans(Z)
      deltanumtor <- deltadentor <- NULL
      deltadentor <- colSums((1-D)*(1-Z))
      deltanumtor <- colSums((1-D)*(1-Z)*dat)
      delta <- deltanumtor/deltadentor
      mu_i_final <- (colSums(log_dat*D)+mis_sum)/(sum(D)+T_i)
      sigma_i_final <- sqrt((colSums(log_dat^2*D)+mis_sum_sq)/(sum(D)+T_i)-mu_i_final^2)
      #end of innter M step
      #inner E (E2) step for Z (ZIP component)
      if (K!=5)
      {
        for (i in 1:dim(dat)[2])
        {
          Z[,i] <- pi_j^(1-D)/(pi_j^(1-D)+(1-pi_j)^(1-D)*exp(-delta[i]*(1-D)))
        }
        Z[dat!=0] <- 0
        #inner E (E2) step for T and Lt (Truncated Normal component augment data)
        beta <- -mu_i_final/sigma_i_final
        T_i <- sum(D)*pnorm(beta)/(1-pnorm(beta))
        mis_sum <- T_i*(mu_i_final-dnorm(-beta)*sigma_i_final/pnorm(beta))
        mis_sum_sq <- T_i*(sigma_i_final^2+mu_i_final^2+sigma_i_final^2*(-beta)*dnorm(-beta)/pnorm(beta)-2*mu_i_final*sigma_i_final*dnorm(beta)/pnorm(beta))
      }
      K <- K+1
    }
    #outer E step
    log_numtor <- apply(matrix(log(truncnorm::dtruncnorm(t(log_dat),a=0,b=Inf, mu_i_final, sigma_i_final)),dim(log_dat)[2],dim(log_dat)[1]), 2, sum)
    log_dentor1 <- dzip(x=dat,delta=delta,pi_j=pi_j)
    #correct for numerical issues
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


#' Produce MIXnorm normalized expression matrix
#'
#' Calls the MIXnorm.mle function to obtain MLE of the mixture model,
#' then produces the normalized expression matrix.
#' @param dat input raw read count matrix. dim of dat = J genes * I samples.
#' @param max_iter maximum number of iterations for the nested EM algorithm default is 20, recommend range (10, 50).
#' @param tol convergency criteria, default is 1e-2, recommend range (1e-5,1).
#' @param appr binary True of False, indicates if the approximate version of normalization should be used.
#' @return A list contains the normalized expression matrix (MIX_normalized_log), proportion of expressed genes (phi) and probabilities of being expressed for all genes (D).
#' @export
MIXnorm <- function(dat,max_iter = 20, tol = 1e-2,appr=T)
{
  MIX_res <- MIXnorm.mle(dat,max_iter=max_iter, tol=tol)
  beta <- (-MIX_res$mu_i_final)/MIX_res$sigma_i_final
  subt <- matrix(rep((MIX_res$mu_i_final+dnorm(beta)*MIX_res$sigma_i_final/(1-pnorm(beta))),each=dim(dat)[1]),dim(dat)[1],dim(dat)[2])
  #MIX_log_normalized returns the normalized data in log scale
  if (!appr)
  {
    MIX_normalized_log <- as.data.frame(t(MIX_res$D*(log(dat + 1)-subt)))
  }
  if (appr)
  {
    MIX_normalized_log <- as.data.frame(t((MIX_res$D>0.5)*(log(dat + 1)-subt)))
  }
  names(MIX_res$D) <- rownames(dat)
  if (MIX_res$exit.code==0)
  {
    cat('Normalization finished successfully \n')
  }
  return(list(MIX_normalized_log = t(MIX_normalized_log), phi = MIX_res$phi, D = MIX_res$D))
}




