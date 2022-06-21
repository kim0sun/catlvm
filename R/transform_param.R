se_logit_fi <- function(scores) {
   score <- do.call(rbind, lapply(scores, do.call, what = rbind))
   fi <- score %*% t(score)
   vcov <- MASS::ginv(fi)
}

get_se <- function(vcov) {
   var <- diag(vcov)
   var[var < 0] = 0
   sqrt(var)
}

se_logit_nlm <- function(logit_par, data, args) {
   lparam <- unlist(logit_par)
   indInf <- is.infinite(lparam) & lparam > 0
   indNegInf <- is.infinite(lparam) & lparam < 0
   nparam <- length(lparam)

   nlm_fit <- nlm(
      floglik, lparam[!(indInf|indNegInf)],
      y = data$y, nobs = args$nobs,
      nvar = args$nvar, ncat = args$ncat,
      nlv = args$nlv, nroot = args$nroot,
      nlink = args$nlink, nlink_unique = args$nlink_unique,
      nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
      root = args$root - 1, ulv = args$u - 1,
      vlv = args$v - 1, cstr_link = args$cstr_link - 1,
      leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
      nclass = args$nclass, nclass_leaf = args$nclass_leaf,
      nclass_u = args$nclass_u, nclass_v = args$nclass_v,
      indInf = indInf, indNegInf = indNegInf, npar = nparam,
      iterlim = 1, hessian = TRUE
   )

   hessian <- nlm_fit$hessian
   vcov <- MASS::ginv(hessian)
   var <- diag(vcov)
   var[var < 0] <- 0

   se = numeric(nparam)
   se[indInf|indNegInf] <- 0
   se[!(indInf|indNegInf)] <- sqrt(var)

   list(vcov = vcov,
        se = splitSE(se, args$ncat, args$nroot,
                     args$nlink_unique, args$nleaf_unique,
                     args$root - 1, args$u - 1, args$v - 1,
                     args$nclass, args$nclass_u,
                     args$nclass_v, args$nclass_leaf))
}

bdiag <- function(x) {
   nr = sum(sapply(x, nrow))
   nc = sum(sapply(x, ncol))
   ans = matrix(0, nr, nc)
   ibegin = 1; iend = 0
   jbegin = 1; jend = 0

   for (m in x){
      iend = iend + nrow(m)
      jend = jend + ncol(m)
      ans[ibegin:iend, jbegin:jend] = m
      ibegin = ibegin + nrow(m)
      jbegin = jbegin + ncol(m)
   }

   ans
}

jac_logistic <- function(x, simplify = TRUE) {
   ex <- exp(x)
   len <- length(ex)

   diag(ex, len, len - 1) - ex ** 2
}

jac_logistic(tau[[1]])
jac_logistic <- function(x) {
   jac_pi <- lapply(pi, jac_logistic)
   jac_tau <- lapply(tau, apply, 2, function(y) list(jac_logistic(y)))
   jac_rho <- lapply(rho, apply, 1, jac_logistic, simplify = 'list')
   apply(tau[[1]], 2, jac_logistic)
}

split_rho <- function(x, k, r) {
   split_k = lapply(rho, matrix, ncol = k)
   lapply(seq(length(split_k)), function(x)
      lapply(split_k[[x]], split, rep(seq(length(ncat[[x]])), ncat[[x]])))
   sweep(split_k, )

   lapply(apply(split_k[[1]], 2, list), lapply, split, )
}

split(1:6, rep(seq(length(ncat[[x]])), ncat[[x]]))
ncat = args$ncat

se_transform <- function(vcov, args) {


   index = rep(1:3, args$npar)
   vcov_pi  <- covmat[index == 1, index == 1]
   vcov_tau <- covmat[index == 2, index == 2]
   vcov_rho <- covmat[index == 3, index == 3]

   pi  <- lapply(log_par$pi, exp)
   tau <- lapply(log_par$tau, exp)
   rho <- lapply(log_par$rho, exp)

   jac_pi <- bdiag(lapply(pi, jacobian))
   jac_tau <- bdiag(lapply(tau, function(x)
      bdiag(apply(x, 2, jacobian, simplify = FALSE))))
   jac_rho <- bdiag(lapply(1:length(rho), function(x)
      bdiag(lapply(split_by(rho[[x]], ncat[[x]]), jacobian))))

   se_pi <- diag(jac_pi %*% vcov_pi %*% t(jac_pi))
   se_pi[se_pi < 0] <- 0
   se_tau <- diag(jac_tau %*% vcov_tau %*% t(jac_tau))
   se_tau[se_tau < 0] <- 0
   se_rho <- diag(jac_rho %*% vcov_rho %*% t(jac_rho))
   se_rho[se_rho < 0] <- 0

   list(pi = sqrt(se_pi),
        tau = sqrt(se_tau),
        rho = sqrt(se_rho))
}
