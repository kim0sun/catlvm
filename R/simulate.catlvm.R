#' Simulate catlvm object
#'
#'
#' @export
simulate.catlvm <- function(
   object, nsim = 500, params = NULL, ...
) {
   args <- object$args
   args$ncat <- if (!is.null(args$ncat)) args$ncat
   else lapply(args$nvar, function(x) rep(2, x))

   if (object$estimated) {
      params <- object$args$log_par
      pi <- params$pi
      tau <- params$tau
      rho <- params$rho
   } else {
      pi <- params$pi
      tau <- params$tau
      rho <- params$rho
      if (length(pi) < args$nroot)
         pi <- lapply(seq_len(args$nroot), function(x) NULL)
      for (r in seq_len(args$nroot)) {
         pi[[r]] <- pi_valid(pi[[r]], args$nclass[args$root[r]], TRUE)
      }
      if (length(tau) < args$nlink_unique)
         tau <- lapply(seq_len(args$nlink_unique), function(x) NULL)
      for (d in seq_len(args$nlink_unique)) {
         tau[[d]] <- tau_valid(tau[[d]], args$nclass_u[d], args$nclass_v[d], TRUE)
      }
      if (length(rho) < args$nleaf_unique)
         rho <- lapply(seq_len(args$nleaf_unique), function(x) NULL)
      for (v in seq_len(args$nleaf_unique)) {
         rho[[v]] <- rho_valid(rho[[v]], args$nclass_leaf[v], args$ncat[[v]], TRUE)
      }
   }

   ysim <- ysim(
      nsim, args$ncat, args$nlv, args$root - 1, args$leaf - 1,
      args$u - 1, args$v - 1, args$cstr_link - 1, args$cstr_leaf - 1,
      args$nroot, args$nleaf, args$nlink, args$nclass,
      pi, tau, rho, TRUE
   )

   # data.name
   y <- data.frame(do.call(cbind, lapply(ysim$y, t)))
   colnames(y) <- unlist(object$model$vars$manifest)

   list(response = y, class = ysim$class, args = args,
        params = output_param(list(pi = pi, tau = tau, rho = rho),
                              object$model, args))
}
