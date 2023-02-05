#' @title Sample Size Calculation for Unstratified Time-to-Event Regression Models
#' @description Calculate the required sample size for testing effects of main variable(s) in unstratified time-to-event (Cox and Fine-Gray) regression models
#'
#' @param alpha Type I error rate
#' @param power Power level
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ Variance for main variable(s); should be a matrix if more than one variable
#' @param R2 Coefficient of determination matrix of main variable(s) given other covariate(s)
#' @param sstype Type of sample size calculation ('N' = patients (default), 'Events' = # events)
#' @param psi Event probability; required for patient-drive sample size (sstype = 'N')
#' @details The sample size in terms of subjects is calculated according to the formula given for GLMs by Theorems 2-3 of Martens, Kim, and Ahn, of the form
#' \deqn{n = \frac{\kappa(\alpha,\pi,p)}{\psi \boldsymbol{\delta}^T \boldsymbol{\Omega}_z (\mathbf{I} - \mathbf{R}^{\otimes 2}) \boldsymbol{\Omega}_z^T \boldsymbol{\delta}},}{n = \kappa(\alpha,\pi,p) / [\psi \delta^T \Omega_z (I - R * R^T) \Omega_z^T \delta],}
#' where \eqn{\delta} is the set of targeted coefficient values for the main variables of interest (log hazard ratios for Cox model and log subdistribution hazard ratios for Fine-Gray model);
#' \eqn{\psi} is the probability of a subject having an event of interest that is not censored;
#' \eqn{\boldsymbol{\Omega}_z} is the Cholesky root of \eqn{Var \mathbf{Z}}{Var Z};
#' \eqn{\mathbf{R}^{\otimes 2}}{R R^T} is the coefficient of determination matrix of the main variables with other covariates;
#' \eqn{\kappa(\alpha,\pi,p)} is the noncentrality parameter for a \eqn{chi_p^2} distribution that solves
#' \deqn{\pi = Prob(\chi_p^2(\kappa) > C_\alpha);}
#' and \eqn{C_\alpha} is the \eqn{1-\alpha} quantile of a central \eqn{chi_p^2} distribution.
#'
#' For an event-driven calculation, the number of events required is calculated according to the formula given for GLMs by Theorem 1 of Martens, Kim, and Ahn, of the form
#' \deqn{e = \frac{\kappa(\alpha,\pi,p)}{\boldsymbol{\delta}^T \boldsymbol{\Omega}_z (\mathbf{I} - \mathbf{R}^{\otimes 2}) \boldsymbol{\Omega}_z^T \boldsymbol{\delta}}.}{n = \kappa(\alpha,\pi,p) / [\delta^T \Omega_z (I - R * R^T) \Omega_z^T \delta].}
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
#' vz = matrix(c(1,0.25,0.25,0.25),nrow=2)
#' rsq = matrix(c(0.5,0.25,0.25,0),nrow=2)
#' pr_event = 0.5
#' ssTTE(level,pow,del,vz,rsq,"Events")     # Event driven SS calculation
#' ssTTE(level,pow,del,vz,rsq,"N",pr_event) # Patient driven SS calculation
ssTTE = function(alpha,power,delta,varZ,R2,sstype='N',psi=NULL) {
  if(tolower(sstype) == 'n' & is.null(psi)) {
    print("Error: psi needs to be specified when sstype='N'")
  }
  else if(tolower(sstype) %in% c('n','events') == 0) {
    print("Error: sstype must equal 'N' or 'Events'")
  }
  else {
    p = length(delta)
    sdZ = chol(varZ)
    numer = uniroot(function(x) pchisq(qchisq(1-alpha,p),p,ncp=x,lower.tail=FALSE)-power,c(0,10^6))$root
    psi = ifelse(tolower(sstype)=='n',psi,1)
    denom = psi*t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2) %*% t(sdZ) %*% delta
    return(as.numeric(numer/denom))
  }
}
