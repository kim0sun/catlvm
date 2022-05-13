##' @export
estimate <- function(object, ...) UseMethod("estimate")

##' @export
estimate.catlvm <- function(
   object, data = parent.frame(), regression = NULL,
   method = "hybrid", smoothing = TRUE, control = catlvm.control(), ...)
{
   if (!is.null(object$data)) {
      data <- object$data
      args <- object$args
   }
   else {
      data <- proc_data(data, object$struct)
      args <- update_args(object$args, data)
      object$data <- data
      object$args <- args
   }

   if (!inherits(control, "catlvm.control")) {
      ctrl <- catlvm.control()
      index <- match(names(control), names(ctrl), 0L)
      ctrl[index] <- control[index > 0]
      control <- ctrl
   }

   init.param <- control$init.param
   is.init <- init_validate(init.param, args)
   for (i in names(is.init)[as.logical(is.init)]) {
      init.param[[i]] <- lapply(init.param[[i]], log)
   }

   if (method == "em") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nedge,
         args$nleaf, args$nleaf_unique,
         args$root - 1, args$tree_index - 1,
         args$u - 1, args$v - 1,
         args$leaf - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         is.init, init.param,
         control$em.iterlim, control$em.tol,
         control$verbose, control$per.iter
      )
      log_par <- em_fit$params
      posterior <- em_fit$posterior
      em.convergence <- em_fit$converged
      nlm.convergence <- NA
      if (control$verbose) cat(".. done.")
   } else if (method == "nlm") {
      params <- par_gnr(args$nobs, args$nvar, args$ncat,
                        args$nroot, args$nedge, args$nleaf_unique,
                        args$root - 1, args$u - 1, args$v - 1,
                        args$nclass, args$nclass_leaf,
                        is.init, init.param)
      lparam <- logit_param(params, args)
      if (control$verbose) cat("\nnlm iteration begin.")
      nlm_fit <- nlm(
         floglik, unlist(lparam),
         y = data$y, nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
         nlv = args$nlv, nroot = args$nroot, nedge = args$nedge,
         nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
         root = args$root - 1, ulv = args$u - 1, vlv = args$v - 1,
         leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
         nclass = args$nclass, nclass_leaf = args$nclass_leaf,
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      log_par <-logit2log(nlm_fit$estimate, args$nobs, args$ncat,
                          args$nroot, args$nedge, args$nleaf_unique,
                          args$root - 1, args$u - 1, args$v - 1,
                          args$nclass, args$nclass_leaf)
      em.convergence <- NA
      nlm.convergence <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.")
   } else if (method == "hybrid") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nedge,
         args$nleaf, args$nleaf_unique,
         args$root - 1, args$tree_index - 1,
         args$u - 1, args$v - 1,
         args$leaf - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         is.init, init.param,
         control$em.iterlim, control$em.tol,
         control$verbose, control$per.iter
      )
      lparam <- logit_param(em_fit$params, args)
      posterior <- em_fit$posterior
      em.convergence <- em_fit$converged
      if (control$verbose) cat(".. done. \nnlm iteration begin.")

      nlm_fit <- nlm(
         floglik, unlist(lparam),
         y = data$y, nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
         nlv = args$nlv, nroot = args$nroot, nedge = args$nedge,
         nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
         root = args$root - 1, ulv = args$u - 1, vlv = args$v - 1,
         leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
         nclass = args$nclass, nclass_leaf = args$nclass_leaf,
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      log_par <- logit2log(nlm_fit$estimate, args$nobs, args$ncat,
                           args$nroot, args$nedge, args$nleaf_unique,
                           args$root - 1, args$u - 1, args$v - 1,
                           args$nclass, args$nclass_leaf)
      nlm.convergence <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.")
   }
   # else if (method == "vem") {
   #    fit <- emFit(data$y, args$nobs, args$nvar, args$ncat,
   #                 args$nlv, args$root - 1, args$leaf - 1,
   #                 args$u - 1, args$v - 1,
   #                 args$tree_index - 1, args$cstr_leaf - 1,
   #                 args$nclass, args$nclass_leaf,
   #                 is.init, init.param,
   #                 control$em.iterlim, control$em.tol,
   #                 control$verbose, control$per.iter)
   # } else if (method == "gibbs") {
   #    fit <- emFit(data$y, args$nobs, args$nvar, args$ncat,
   #                 args$nlv, args$root - 1, args$leaf - 1,
   #                 args$u - 1, args$v - 1,
   #                 args$tree_index - 1, args$cstr_leaf - 1,
   #                 args$nclass, args$nclass_leaf,
   #                 is.init, init.param,
   #                 control$em.iterlim, control$em.tol,
   #                 control$verbose, control$per.iter)
   # }
   logit_par <- logit_param(log_par, args)
   se_par <- se_param(logit_par, data, args)
   posterior <- calcPost(
      log_par, data$y, args$nobs, args$nvar, args$ncat,
      args$nlv, args$nroot, args$nedge, args$nleaf, args$nleaf_unique,
      args$root - 1, args$tree_index - 1, args$u - 1, args$v - 1,
      args$leaf - 1, args$cstr_leaf - 1, args$nclass, args$nclass_leaf
   )

   fit = list()
   fit$estimates = list(
      par = output_param(log_par, object$struct, args),
      logit = output_logit(logit_par, object$struct, args),
      se_logit = output_se(se_par, object$struct, args)
   )
   fit$loglik <- floglik(
      unlist(logit_par), data$y, args$nobs, args$nvar, args$ncat,
      args$nlv, args$nroot, args$nedge, args$nleaf, args$nleaf_unique,
      args$root - 1, args$u - 1, args$v - 1,
      args$leaf - 1, args$cstr_leaf - 1,
      args$nclass, args$nclass_leaf
   )
   fit$posterior <- output_posterior(posterior, object$struct, data)
   fit$convergence <- c(EM = em.convergence, nlm = nlm.convergence)
   object$fit <- fit

   class(object) <- "catlvm.fit"
   object
}
