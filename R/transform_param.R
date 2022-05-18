logit_param <- function(param, args) {
   pi <- list()
   tau <- list()
   rho <- list()

   for (r in seq_len(args$nroot)) {
      nclass <- args$nclass[args$root[r]]
      pi[[r]] <- logit_pi(param$pi[[r]], nclass)
   }

   for (d in seq_len(args$nedge)) {
      nk <- args$nclass[args$u[d]]
      nl <- args$nclass[args$v[d]]
      tau[[d]] <- logit_tau(param$tau[[d]], nk, nl)
   }

   for (v in seq_len(args$nleaf_unique)) {
      nclass <- args$nclass_leaf[v]
      rho[[v]] <- logit_rho(param$rho[[v]], args$ncat[[v]], nclass)
   }

   return(list(pi = pi, tau = tau, rho = rho))
}

se_logit_par <- function(logit_par, data, args) {
   nlm_fit <- nlm(
      floglik, unlist(logit_par),
      y = data$y, nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
      nlv = args$nlv, nroot = args$nroot, nedge = args$nedge,
      nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
      root = args$root - 1, ulv = args$u - 1, vlv = args$v - 1,
      leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
      nclass = args$nclass, nclass_leaf = args$nclass_leaf,
      iterlim = 1, hessian = TRUE
   )

   hessian <- nlm_fit$hessian
   vcov <- MASS::ginv(hessian)
   var <- diag(vcov)
   var[var < 0] <- 0
   se <- sqrt(var)

   list(vcov = vcov,
        se = splitSE(se, args$ncat, args$nroot,
                     args$nedge, args$nleaf_unique,
                     args$root - 1, args$u - 1, args$v - 1,
                     args$nclass, args$nclass_leaf))
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

jacobian <- function(x) {
   len <- length(x)
   diag(x, len + 1, len) - x ** 2
}

split_by <- function(x, y) split(x, rep(1:length(y), y))

se_param <- function(covmat, log_par, args) {
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
