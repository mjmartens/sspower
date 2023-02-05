#' @title Sample Size Calculation for Stratified Time-to-Event Regression Models
#' @description Calculate the required sample size for testing effects of main variable(s) in stratified time-to-event (Cox and Fine-Gray) regression models
#'
#' @param alpha Type I error rate
#' @param power Power level
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ 3-dimensional array of variance matrices for main variable(s) in strata; each slice corresponds to a stratum
#' @param R2 Array of coefficient of determination matrices of main variable(s) given other covariate(s) for strata; each slice corresponds to a stratum
#' @param psi Vector of event probabilities for strata
#' @param rho Vector of strata probabilities
#' @param sstype Type of sample size calculation ('N' = patients (default), 'Events' = # events)
#' @details The sample size in terms of subjects is calculated according to the formula given for GLMs by Theorems 2-3 of Martens, Kim, and Ahn, of the form
#' \deqn{n = \frac{\kappa(\alpha,\pi,p)}{\sum_{k=1}^K \rho_k \psi_k \boldsymbol{\delta}^T \boldsymbol{\Omega}_{kz} (\mathbf{I} - \mathbf{R}_k^{\otimes 2}) \boldsymbol{\Omega}_{kz}^T \boldsymbol{\delta}},}{n = \kappa(\alpha,\pi,p) / [\sum_{k=1}^K \rho_k \psi_k \delta^T \Omega_kz (I - R_k * R_k^T) \Omega_kz^T \delta],}
#' where \eqn{\delta} is the set of targeted coefficient values for the main variables of interest (log hazard ratios for Cox model and log subdistribution hazard ratios for Fine-Gray model);
#' \eqn{\psi}_k is the probability of a subject in stratum k having an event of interest that is not censored;
#' \eqn{\rho}_k is the probability of being in stratum k;
#' \eqn{\boldsymbol{\Omega}_{kz}} is the Cholesky root of \eqn{Var \mathbf{Z}}{Var Z} in stratum k;
#' \eqn{\mathbf{R}_k^{\otimes 2}}{R_k R_k^T} is the coefficient of determination matrix of the main variables with other covariates in stratum k;
#' \eqn{\kappa(\alpha,\pi,p)} is the noncentrality parameter for a\eqn{chi_p^2} distribution that solves
#' \deqn{\pi = Prob(\chi_p^2(\kappa) > C_\alpha);}
#' and \eqn{C_\alpha} is the \eqn{1-\alpha} quantile of a central \eqn{chi_p^2} distribution.
#'
#' For an event-driven calculation, the number of events required is calculated according to the formula given for GLMs by Theorem 1 of Martens, Kim, and Ahn, of the form
#' \deqn{e = \frac{\kappa(\alpha,\pi,p) \sum_{k=1}^K \rho_k \psi_k}{\sum_{k=1}^K \rho_k \psi_k \boldsymbol{\delta}^T \boldsymbol{\Omega}_{kz} (\mathbf{I} - \mathbf{R}_k^{\otimes 2}) \boldsymbol{\Omega}_{kz}^T \boldsymbol{\delta}},}{n = \kappa(\alpha,\pi,p) \sum_{k=1}^K \rho_k \psi_k / [\sum_{k=1}^K \rho_k \psi_k \delta^T \Omega_kz (I - R_k * R_k^T) \Omega_kz^T \delta],}
#' @return Required sample size to detect targeted effect sizes given specified type I error rate and power
#' @export
#'
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#' @examples
#' ## Example of calculating sample size to attain 80\% power for a 5\% significance level test
#' level = 0.05
#' pow = 0.80
#' del = c(0.4,0.3)
#' vz = array(0,c(2,2,2))
#' vz[,,1] = diag(rep(0.25,2))
#' vz[,,2] = diag(rep(1,2))
#' pr_event = c(0.2,0.3)
#' pr_strat = c(0.3,0.7)
#' rsq = array(0,c(2,2,2))
#' rsq[,,1] = diag(rep(0,2))
#' rsq[,,2] = matrix(c(0.5,0.25,0.25,0),nrow=2)
#' ssTTEstrat(level,pow,del,vz,rsq,pr_event,pr_strat,"Events")  # Event driven SS calculation
#' ssTTEstrat(level,pow,del,vz,rsq,pr_event,pr_strat,"N")       # Patient driven SS calculation
ssTTEstrat = function(alpha,power,delta,varZ,R2,psi,rho,sstype='N') {
  if(tolower(sstype) %in% c('n','events') == 0) {
    print("Error: sstype must equal 'N' or 'Events'")
  }
  else {
    p = length(delta)
    numer = uniroot(function(x) pchisq(qchisq(1-alpha,p),p,ncp=x,lower.tail=FALSE)-power,c(0,10^6))$root
    numer = numer*sum(rho*psi)^(tolower(sstype)=='events')
    denom = 0
    for(k in 1:dim(varZ)[3]){
      sdZ = chol(varZ[,,k])
      denom = denom + psi[k] * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2[,,k]) %*% t(sdZ) %*% delta
    }
    return(as.numeric(numer/denom))
  }
}
