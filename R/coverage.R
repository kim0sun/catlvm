#' @export
coverage <- function(
   object, niter = 100, parm, nsim = 500,
   level = 0.95, verbose = TRUE
) {
   if (object$estimated)
      parm <- object$estimates$param
   else if (!missing(parm)) {
      sim <- object %>% simulate(nsim = nsim, param = parm)
      parm <- sim$params
   } else {
      sim <- object %>% simulate(nsim = nsim)
      parm <- sim$params
   }

   cover <- par <- unlist(parm)
   cover[] <- 0

   for (i in 1:niter) {
      if (verbose && i %% 50 == 0) cat(".")
      sim <- object %>% simulate(nsim = nsim, params = parm)
      fit <- object %>% estimate(
         data = sim$response,
         control = list(verbose = FALSE, init.param = parm)
      )
      ci <- confint.catlvm(fit, level = level)
      cover <- cover + as.numeric(par > ci[, 1] & par < ci[, 2])
   }
   if (verbose) cat("DONE. \n")
   cover / niter
}
