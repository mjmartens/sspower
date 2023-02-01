#' @title Power Calculation for Unstratified Time-to-Event Regression Models
#' @description Calculate the power attained for testing effects of main variable(s) in unstratified time-to-event (Cox and Fine-Gray) regression models
#'
#' @param n Number of subjects
#' @param alpha Type I error rate
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ Variance for main variable(s); should be a matrix if more than one variable
#' @param R2 Coefficient of determination matrix of main variable(s) given other covariate(s)
#' @param psi Event probability
#' @return Power to detect targeted effect sizes given specified sample size and type I error rate
#' @export
#'
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#' @examples
#' ## Example of calculating sample size to attain 80\% power for a 5\% significance level test
#' num = 200
#' level = 0.05
#' del = c(0.4,0.3)
#' vz = matrix(c(1,0.25,0.25,0.25),nrow=2)
#' rsq = matrix(c(0.5,0.25,0.25,0),nrow=2)
#' pr_event = 0.5
#' powerTTE(num,level,del,vz,rsq,pr_event)
powerTTE = function(n,alpha,delta,varZ,R2,psi) {
  p = length(delta)
  sdZ = chol(varZ)
  kappa = n*psi * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2) %*% t(sdZ) %*% delta
  val = pchisq(qchisq(1-alpha,p),p,ncp=kappa,lower.tail=FALSE)
  return(val)
}
