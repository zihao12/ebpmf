#' @title Empirical Bayes Poisson Matrix Factorization
#' @description Uses Empirical Bayes to fit the model \deqn{X_{ij}  ~ Poi(\sum_k L_{ik} F_{jk})} with \deqn{L_{.k} ~ g_k()}
#' with g_k being either Mixture of Exponential, or Point Gamma
#' @import mixsqp
#' @import ebpm

#' @details The model is fit in 2 stages: i) estimate \eqn{g} by maximum likelihood (over pi_k)
#' ii) Compute posterior distributions for \eqn{\lambda_j} given \eqn{x_j,\hat{g}}.
#' @param X count matrix (dim(X) = c(n, p)).
#' @param K number of topics
#' @param m multiplicative parameter for selecting grid in "ebpm::ebpm_exponential_mixture"

#' @param  maxiter.out  maximum iterations in the outer  loop
#' @param  maxiter.int  maximum iterations in the inner loop
#' @param seed random seed
#'
#'
#' @return A list containing elements:
#'     \describe{
#'       \item{\code{qls_mean}}{A n by k matrix: Approximate posterior mean for L}
#'       \item{\code{qls_mean_log}}{A n by k matrix: Approximate posterior log mean for L}
#'       \item{\code{gls}}{A list of K elements, each element is the estimated prior for the kth column of L}
#'       \item{\code{qfs_mean}}{A p by k matrix: Approximate posterior mean for F}
#'       \item{\code{qfs_mean_log}}{A p by k matrix: Approximate posterior log mean for F}
#'       \item{\code{gfs}}{A list of K elements, each element is the estimated prior for the kth column of F}
#'      }
#'
#'
#' @examples
#' To add

#'
#' @export


ebpmf_exponential_mixture <- function(X, K, m = 2, maxiter.out = 10, maxiter.int = 1, seed = 123){
  set.seed(123)
  start = proc.time()
  qg = initialize_qg(X, K)
  runtime_init = (proc.time() - start)[[3]]
  runtime_rank1 = 0
  runtime_ez = 0

  for(iter in 1:maxiter.out){
    #print(sprintf("iter: %d", iter))
    for(k in 1:K){
      ## get row & column sum of <Z_ijk>
      start = proc.time()
      Ez = get_Ez(X, qg, k, K)
      runtime_ez = runtime_ez + (proc.time() - start)[[3]]
      ## update q, g
      start = proc.time()
      tmp = ebpmf_rank1_exponential_helper(Ez$rsum,Ez$csum,NULL,m, maxiter.int)
      runtime_rank1 = runtime_rank1 + (proc.time() - start)[[3]]
      qg = update_qg(tmp, qg, k)
    }
  }
  print("summary of  runtime:")
  print(sprintf("init           : %f", runtime_init))
  print(sprintf("Ez     per time: %f", runtime_ez/(iter*K)))
  print(sprintf("rank1  per time: %f", runtime_rank1/(iter*K)))
  return(qg)
}

## ================== helper functions ==================================
## for each pair of l, f, give them 1/k of the row & col sum
initialize_qg <- function(X, K, seed = 123){
  n = nrow(X)
  p = ncol(X)
  set.seed(seed)
  X_rsum = rowSums(X)
  X_csum = colSums(X)
  prob_r = replicate(n, rdirichlet(1,replicate(K, 1/K)))[1,,] ## K by n
  prob_c = replicate(p, rdirichlet(1,replicate(K, 1/K)))[1,,] ## K  by p
  rsums = matrix(replicate(K*n,0), nrow = K)
  csums = matrix(replicate(K*p,0), nrow = K)
  for(i in  1:n){
    if(X_rsum[i] == 0){rsums[,i] = replicate(K, 0)}
    else{rsums[,i] = rmultinom(1, X_rsum[i],prob_r[,i])}
  }
  for(j in  1:p){
    if(X_csum[j] == 0){csums[,j] = replicate(K, 0)}
    else{csums[,j] = rmultinom(1, X_csum[j],prob_c[,j])}
  }
  qg = list(qls_mean = matrix(replicate(n*K, 0), ncol =  K), qls_mean_log = matrix(replicate(n*K, 0), ncol =  K), gls = replicate(K, list(NaN)),
            qfs_mean = matrix(replicate(p*K, 0), ncol =  K), qfs_mean_log = matrix(replicate(p*K, 0), ncol =  K), gfs = replicate(K, list(NaN))
  )
  for(k in 1:K){
    qg_ = ebpmf_rank1_exponential_helper(rsums[k,], csums[k, ], init = NULL, m = 2, maxiter = 1)
    qg   = update_qg(qg_, qg, k)
  }
  return(qg)
}

## compute the row & col sum of <Z_ijk> for a given k
get_Ez <- function(X, qg, k, K){
  n = nrow(X)
  p = ncol(X)
  psi = array(dim = c(n, p, K))
  ## get <ln l_ik> + <ln f_jk>
  for(d in 1:K){
    psi[,,d] = outer(qg$qls_mean_log[,d], qg$qfs_mean_log[,d], "+")
  }
  ## do softmax
  #browser()
  psi = softmax3d(psi)
  Ez = as.vector(psi)*as.vector(X)
  dim(Ez) = dim(psi)
  return(list(rsum = rowSums(Ez[,,k]), csum = colSums(Ez[,,k])))
}

softmax3d <- function(x){
  score.exp <- exp(x)
  probs <-as.vector(score.exp)/as.vector(rowSums(score.exp,dims=2))
  dim(probs) <- dim(x)
  return(probs)
}

update_qg <- function(tmp, qg, k){
  qg$qls_mean[,k] = tmp$ql$mean
  qg$qls_mean_log[,k] = tmp$ql$mean_log
  qg$qfs_mean[,k] = tmp$qf$mean
  qg$qfs_mean_log[,k] = tmp$qf$mean_log
  qg$gls[[k]] = tmp$gl
  qg$gfs[[k]] = tmp$gf
  return(qg)
}

ebpmf_rank1_exponential_helper <- function(X_rowsum,X_colsum, init = NULL, m = 2, maxiter = 1){
  if(is.null(init)){init = list(mean = runif(length(X_rowsum), 0, 1))}
  ql = init
  for(i in 1:maxiter){
    ## update q(f), g(f)
    sum_El = sum(ql$mean)
    tmp = ebpm::ebpm_exponential_mixture(x = X_colsum, s = replicate(p,sum_El), m = m)
    qf = tmp$posterior
    gf = tmp$fitted_g
    ll_f = tmp$log_likelihood
    ## update q(l), g(l)
    sum_Ef = sum(qf$mean)
    tmp = ebpm_exponential_mixture(x = X_rowsum, s = replicate(n,sum_Ef), m = m)
    ql = tmp$posterior
    gl = tmp$fitted_g
    ll_l = tmp$log_likelihood
    qg = list(ql = ql, gl = gl, qf = qf, gf = gf, ll_f = ll_f, ll_l = ll_l)
  }
  return(qg)
}















