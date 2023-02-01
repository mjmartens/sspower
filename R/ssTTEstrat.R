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
#' @return Required sample size to detect targeted effect sizes with specified type I error rate and power
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
