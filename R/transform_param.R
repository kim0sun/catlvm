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

se_param <- function(logit_par, data, args) {
   nlm_fit <- nlm(
      llFit, unlist(logit_par),
      y = data$y, nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
      nlv = args$nlv, nroot = args$nroot, nedge = args$nedge,
      nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
      root = args$root - 1, ulv = args$u - 1, vlv = args$v - 1,
      leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
      nclass = args$nclass, nclass_leaf = args$nclass_leaf,
      iterlim = 1, hessian = TRUE
   )

   hessian <- nlm_fit$hessian
   var <- diag(MASS::ginv(hessian))
   var[var < 0] <- 0
   se <- sqrt(var)

   splitSE(se, args$ncat, args$nroot, args$nedge, args$nleaf_unique,
           args$root - 1, args$u - 1, args$v - 1,
           args$nclass, args$nclass_leaf)
}
