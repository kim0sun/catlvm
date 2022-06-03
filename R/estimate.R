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
      data <- proc_data(data, object$model)
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
   if (is.null(init.param)) init.param <- list()
   is.init <- init_validate(init.param, args)
   for (i in names(is.init)[as.logical(is.init)]) {
      init.param[[i]] <- lapply(init.param[[i]], log)
   }

   if (method == "em") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$tree_index - 1, args$root - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         is.init, init.param,
         control$em.iterlim, control$em.tol,
         control$verbose, control$per.iter
      )
      log_par <- em_fit$params
      posterior <- em_fit$posterior
      em.convergence <- em_fit$converged
      nlm.convergence <- NA
      if (control$verbose) cat(".. done.\n")
   } else if (method == "nlm") {
      params <- par_gnr(args$nobs, args$nvar, args$ncat,
                        args$nroot, args$nlink_unique, args$nleaf_unique,
                        args$root - 1, args$u - 1, args$v - 1,
                        args$nclass, args$nclass_leaf,
                        args$nclass_u, args$nclass_v,
                        is.init, init.param)
      lparam <- unlist(logit_param(params, args))
      nparam <- length(lparam)
      indInf <- rep(FALSE, nparam)
      indNegInf <- rep(FALSE, nparam)
      if (control$verbose) cat("nlm iteration begin.\n")
      nlm_fit <- nlm(
         floglik, lparam, y = data$y, nobs = args$nobs,
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
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      log_par <-logit2log(nlm_fit$estimate, args$nobs, args$ncat,
                          args$nroot, args$nlink_unique, args$nleaf_unique,
                          args$root - 1, args$u - 1, args$v - 1,
                          args$nclass, args$nclass_leaf,
                          args$nclass_u, args$nclass_v)
      em.convergence <- NA
      nlm.convergence <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.\n")
   } else if (method == "hybrid") {
      if (control$verbose) cat("EM iteration begin.\n")
      em_fit <- emFit(
         data$y, args$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$tree_index - 1, args$root - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v,
         is.init, init.param,
         control$em.iterlim, control$em.tol,
         control$verbose, control$per.iter
      )
      lparam <- unlist(logit_param(em_fit$params, args))
      posterior <- em_fit$posterior
      em.convergence <- em_fit$converged
      if (control$verbose) cat(".. done. \nnlm iteration begin.\n")

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
         iterlim = control$nlm.iterlim, steptol = control$nlm.tol
      )
      logit_par <- numeric(nparam)
      logit_par[!(indInf|indNegInf)] <- nlm_fit$estimate
      logit_par[indInf] <- Inf
      logit_par[indNegInf] <- -Inf

      log_par <- logit2log(logit_par, args$nobs, args$ncat,
                           args$nroot, args$nlink_unique, args$nleaf_unique,
                           args$root - 1, args$u - 1, args$v - 1,
                           args$nclass, args$nclass_leaf,
                           args$nclass_u, args$nclass_v)
      nlm.convergence <- nlm_fit$code < 3
      if (control$verbose) cat(".. done.\n")
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

   etc <- calcModel(
      log_par, data$y, args$nobs, args$nvar, args$ncat,
      args$nlv, args$nroot, args$nlink, args$nleaf,
      args$nlink_unique, args$nleaf_unique,
      args$root - 1, args$tree_index - 1,
      args$u - 1, args$v - 1, args$leaf - 1,
      args$cstr_link - 1, args$cstr_leaf - 1,
      args$nclass, args$nclass_leaf,
      args$nclass_u, args$nclass_v
   )

   logit_par <- logit_param(log_par, args)
   se_logit <- se_logit_par(logit_par, data, args)

   object$estimates = list(
      param = output_param(log_par, object$model, args),
      logit = output_logit(logit_par, object$model, args),
      se_par = NULL,
      se_logit = output_logit(se_logit$se, object$model, args)
   )
   object$llik <- etc$ll
   object$lambda <- etc$lambda[args$root]
   object$posterior <- output_posterior(etc$post, object$model, data)
   object$convergence <- c(EM = em.convergence, nlm = nlm.convergence)
   object$estimated <- TRUE

   class(object) <- "catlvm"
   object
}


