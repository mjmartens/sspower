#' @title Sample Size Calculation for Generalized Linear Models
#'
#' @description Calculate the required sample size for testing effects of main variable(s) in a generalized linear model (GLM)
#'
#' @param alpha Type I error rate
#' @param power Power level
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ Variance for main variable(s); should be a matrix if more than one variable
#' @param R2 Coefficient of determination matrix of main variable(s) given other covariate(s)
#' @param f00 Expected residual variance under null, E[Var(Y|X)|H0]
#'
#' @details The GLM considered has the form
#' \deqn{h\{E(Y|\mathbf{Z},\mathbf{X})\} = \beta_0 + \boldsymbol{\beta}^T \mathbf{Z} + \boldsymbol{\gamma}^T \mathbf{X},}{h\{E(Y|Z,X)\} = \beta_0 + \beta^T Z + \gamma^T X,}
#' where h is the link function, \eqn{\mathbf{Z}}{Z} is the main variable(s) of interest, \eqn{\mathbf{X}}{X} is the other covariate(s), \eqn{\beta_0} is the intercept term, \eqn{\beta} is the set of coefficients for \eqn{\mathbf{Z}}{Z}, and \eqn{\gamma} is the set of coefficients for \eqn{\mathbf{X}}{X}.
#'
#' The sample size is calculated according to the formula given by Theorem 1 of Martens, Kim, and Ahn, of the form
#' \deqn{n = \frac{\kappa(\alpha,\pi,p)}{f_{00} \boldsymbol{\delta}^T \boldsymbol{\Omega}_z (\mathbf{I} - \mathbf{R}^{\otimes 2}) \boldsymbol{\Omega}_z^T \boldsymbol{\delta}},}{n = \kappa(\alpha,\pi,p) / [f_00 \delta^T \Omega_z (I - R * R^T) \Omega_z^T \delta],}
#' where \eqn{\pi} is the power for the test,
#' \eqn{f_{00}}{f_00} is the residual variance \eqn{E[Var(Y|X)|H_0]} under the null hypothesis \eqn{H_0: \beta = 0},
#' \eqn{\delta} is the set of targeted coefficient values for \eqn{\mathbf{Z}}{Z},
#' \eqn{\boldsymbol{\Omega}_z}{\Omega_z} is the Cholesky root of \eqn{Var \mathbf{Z}}{Var Z},
#' \eqn{\mathbf{R}^{\otimes 2}}{R R^T} is the coefficient of determination matrix of \eqn{\mathbf{Z}}{Z} with \eqn{\mathbf{X}}{X},
#' \eqn{\kappa(\alpha,\pi,p)} is the noncentrality parameter for a \eqn{chi_p^2} distribution that solves
#' \deqn{\pi = Prob(\chi_p^2(\kappa) > C_\alpha),}
#' and \eqn{C_\alpha} is the \eqn{1-\alpha} quantile of a central \eqn{chi_p^2} distribution.
#'
#' @return Required sample size to detect targeted effect sizes given specified type I error rate and power
#'
#' @export
#'
#' @references Nelder, J. A. and Wedderburn, R. W. (1972). Generalized linear models. \emph{Journal of the Royal Statistical Society: Series A (General)} \strong{135}, 370–384.
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#'
#' @examples
#' ## Example of calculating sample size to attain 80\% power for a 5\% significance level test
#' level = 0.05
#' pow = 0.80
#' del = c(0.5,0.3)
#' vz = diag(c(0.5,0.25))
#' rsq = diag(c(0,0.5))
#' f0 = 0.20
#' ssGLM(level,pow,del,vz,rsq,f0)
ssGLM = function(alpha,power,delta,varZ,R2,f00) {
  p = length(delta)
  sdZ = chol(varZ)
  numer = uniroot(function(x) pchisq(qchisq(1-alpha,p),p,ncp=x,lower.tail=FALSE)-power,c(0,10^6))$root
  denom = f00 * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2) %*% t(sdZ) %*% delta
  return(as.numeric(numer/denom))
}
