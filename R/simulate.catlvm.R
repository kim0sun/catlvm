##' @export
simulate.catlvm <- function(
   object, nsim = 500,
   params = NULL, ...
) {
   nlv <- object$args$nlv
   root <- object$args$root
   leaf <- object$args$leaf
   ulv <- object$args$u
   vlv <- object$args$v

   cstr_root <- object$args$cstr_root
   cstr_leaf <- object$args$cstr_leaf

   nclass <- object$args$nclass
   nclass_leaf <- object$args$nclass_leaf

   nroot <- object$args$nroot
   nleaf <- object$args$nleaf
   nedge <- object$args$nedge
   nleaf_unique <- object$args$nleaf_unique

   nvar <- object$args$nvar
   ncat <- if (!is.null(object$args$ncat)) obejct$args$ncat
   else lapply(nvar, function(x) rep(2, x))

   pi <- params$pi
   tau <- params$tau
   rho <- params$rho
   if (length(pi) < nroot)
      pi <- lapply(seq_len(nroot), function(x) NULL)
   for (r in seq_len(nroot)) {
      pi[[r]] <- pi_valid(pi[[r]], nclass[root[r]], TRUE)
   }
   if (length(tau) < nedge)
      tau <- lapply(seq_len(nedge), function(x) NULL)
   for (d in seq_len(nedge)) {
      tau[[d]] <- tau_valid(tau[[d]], nclass[ulv[d]], nclass[vlv[d]], TRUE)
   }
   if (length(rho) < nleaf_unique)
      rho <- lapply(seq_len(nleaf_unique), function(x) NULL)
   for (v in seq_len(nleaf_unique)) {
      rho[[v]] <- rho_valid(rho[[v]], nclass_leaf[v], ncat[[v]], TRUE)
   }

   ysim <- ysim(nsim, nlv, root - 1, leaf - 1, ulv - 1, vlv - 1,
                nclass, nroot, nleaf, nedge, ncat, cstr_leaf - 1,
                pi, tau, rho, TRUE)

   # data.name
   y <- data.frame(do.call(cbind, lapply(ysim$y, t)))
   colnames(y) <- unlist(object$struct$vars$manifest)

   args <- object$args

   list(response = y, class = ysim$class, args = args,
        params = list(pi = lapply(pi, exp),
                      tau = lapply(tau, exp),
                      rho = lapply(rho, exp)))
}
