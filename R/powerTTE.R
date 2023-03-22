#' @title Power Calculation for Unstratified Time-to-Event Regression Models
#'
#' @description Calculate the power attained for testing effects of main variable(s) in unstratified time-to-event (Cox and Fine-Gray) regression models
#'
#' @param n Number of subjects
#' @param alpha Type I error rate
#' @param delta Vector of targeted effect size(s) for main variable(s). These correspond to values of the regression coefficients you want to detect.
#' @param varZ Variance for main variable(s); should be a matrix if more than one variable
#' @param R2 Coefficient of determination matrix of main variable(s) given other covariate(s)
#' @param psi Event probability

#' @details The Cox and Fine-Gray models considered have the form
#' \deqn{\lambda(t|\mathbf{Z},\mathbf{X}) = \lambda_0(t) \exp(\boldsymbol{\beta}^T \mathbf{Z} + \boldsymbol{\gamma}^T \mathbf{X}),}{\lambda(t|Z,X) = \lambda_0(t) \exp(\beta^T Z + \gamma^T X),}
#' where \eqn{\lambda} is the hazard function for the Cox model and subdistribution hazard function for the Fine-Gray model,
#' \eqn{\lambda_0} is the baseline hazard/subdistribution hazard function,
#' \eqn{\mathbf{Z}}{Z} is the main variable(s) of interest, \eqn{\mathbf{X}}{X} is the other covariate(s), \eqn{\beta} is the set of coefficients for \eqn{\mathbf{Z}}{Z}, and \eqn{\gamma} is the set of coefficients for \eqn{\mathbf{X}}{X}.
#'
#' The power is calculated according to the formula given by Theorems 2-3 of Martens, Kim, and Ahn, of the form
#' \deqn{power = Prob(\chi_p^2(\kappa) > C_\alpha),}
#' where \deqn{\kappa = n \psi \boldsymbol{\delta}^T \boldsymbol{\Omega}_z (\mathbf{I} - \mathbf{R}^{\otimes 2}) \boldsymbol{\Omega}_z^T \boldsymbol{\delta},}{\kappa = n \psi \delta^T \Omega_z (I - R * R^T) \Omega_z^T \delta,}
#' \eqn{\psi} is the probability of a subject having an event of interest that is not censored,
#' \eqn{\delta} is the set of targeted GLM coefficient values for \eqn{\mathbf{Z}}{Z},
#' \eqn{\boldsymbol{\Omega}_z}{\Omega_z} is the Cholesky root of \eqn{Var \mathbf{Z}}{Var Z},
#' \eqn{\mathbf{R}^{\otimes 2}}{R R^T} is the coefficient of determination matrix of \eqn{\mathbf{Z}}{Z} with \eqn{\mathbf{X}}{X},
#' and \eqn{C_\alpha} is the \eqn{1-\alpha} quantile of a central \eqn{chi_p^2} distribution.
#'
#' @return Power to detect targeted effect sizes given specified sample size and type I error rate
#'
#' @export
#'
#' @references Cox, D. R. (1972). Regression and life tables (with discussion). \emph{Journal of the Royal Statistical Society Series B} \strong{34}, 187–220.
#' @references Fine, J. P. and Gray, R. J. (1999). A proportional hazards model for the subdistribution of a competing risk. \emph{Journal of the American Statistical Association} \strong{94}, 496–509.
#' @references Martens, M.J. and Logan, B.L. (2020). A Unified Approach to Sample Size and Power Determination for Testing Parameters in Generalized Linear and Time to Event Regression Models. \emph{Statistics in Medicine} \strong{40(5)}, 1121-1132.
#' @references Martens, M.J., Kim, S., and Ahn, K.W. (2023+). Sample Size & Power Determination for Evaluating Multiple Parameters in Nonlinear Regression Models with Potential Stratification. \emph{Biometrics}, draft under review.
#'
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
  sdZ = t(chol(varZ))
  kappa = n*psi * t(delta) %*% sdZ %*% (diag(rep(1,p)) - R2) %*% t(sdZ) %*% delta
  val = pchisq(qchisq(1-alpha,p),p,ncp=kappa,lower.tail=FALSE)
  return(val)
}
