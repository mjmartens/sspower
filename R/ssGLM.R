#' @title Sample Size Calculation for Generalized Linear Models
#' @description Calculate the required sample size for testing main variable(s) in a generalized linear model (GLM)
#'
#' @param alpha Type I error rate
#' @param power Power level
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ Variance for main variable(s); should be a matrix if more than one variable
#' @param f00 Expected residual variance under null, E[Var(Y|X)|H0]
#' @param R2 Coefficient of determination matrix of main variable(s) given other covariate(s)
#'
#' @return Required sample size to detect targeted effect sizes with specified type I error rate and power
#' @export
#'
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#' @examples
#' ## Example of calculating sample size to attain 80\% power for a 5\% significance level test
#' level = 0.05
#' pow = 0.80
#' del = c(0.5,0.3)
#' f0 = 0.20
#' vz = diag(c(0.5,0.25))
#' rsq = diag(c(0,0.5))
#' ssGLM(level,pow,del,vz,f0,rsq)
ssGLM = function(alpha,power,delta,varZ,f00,R2) {
  p = length(delta)
  sdZ = chol(varZ)
  numer = uniroot(function(x) pchisq(qchisq(1-alpha,p),p,ncp=x,lower.tail=FALSE)-power,c(0,10^6))$root
  denom = f00 * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2) %*% t(sdZ) %*% delta
  return(as.numeric(numer/denom))
}
