#' @title Power Calculation for Stratified Time-to-Event Regression Models
#' @description Calculate the power attained for testing effects of main variable(s) in stratified time-to-event (Cox and Fine-Gray) regression models
#'
#' @param n Number of subjects
#' @param alpha Type I error rate
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ 3-dimensional array of variance matrices for main variable(s) in strata; each slice corresponds to a stratum
#' @param R2 Array of coefficient of determination matrices of main variable(s) given other covariate(s) for strata; each slice corresponds to a stratum
#' @param psi Vector of event probabilities for strata
#' @param rho Vector of strata probabilities
#' @return Power to detect targeted effect sizes given specified sample size and type I error rate
#' @export
#'
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#' @examples
#' ## Example of calculating sample size to attain 80\% power for a 5\% significance level test
#' num = 200
#' level = 0.05
#' del = rep(0.4,2)
#' vz = array(0,c(2,2,2))
#' vz[,,1] = diag(rep(0.25,2))
#' vz[,,2] = diag(rep(1,2))
#' pr_event = c(0.4,0.3)
#' pr_strat = c(0.3,0.7)
#' rsq = array(0,c(2,2,2))
#' rsq[,,1] = diag(rep(0,2))
#' rsq[,,2] = matrix(c(0.5,0.25,0.25,0),nrow=2)
#' powerTTEstrat(num,level,del,vz,rsq,pr_event,pr_strat)
powerTTEstrat = function(n,alpha,delta,varZ,R2,psi,rho) {
  p = length(delta)
  kappa = 0
  for(k in 1:dim(varZ)[3]){
    sdZ = chol(varZ[,,k])
    kappa = kappa + psi[k] * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2[,,k]) %*% t(sdZ) %*% delta
  }
  kappa = n*kappa
  val = pchisq(qchisq(1-alpha,p),p,ncp=kappa,lower.tail=FALSE)
  return(val)
}