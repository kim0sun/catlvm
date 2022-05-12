##' @export
estimate <- function(object, ...) UseMethod("estimate")

##' @export
estimate.catlvm <- function(
   object, data = parent.frame(), regression = NULL,
   method = "em", smoothing = TRUE, control = catlvm.control(), ...)
{
   if (!is.null(object$data)) {
      data <- object$data
      args <- object$args
   }
   else {
      data <- proc_data(data, object$struct)
      args <- args_return(object$struct, data)
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
      fit <- emFit(data$y, args$nobs, args$nvar, args$ncat,
                   args$nlv, args$root - 1, args$leaf - 1,
                   args$u - 1, args$v - 1,
                   args$tree_index - 1, args$cstr_leaf - 1,
                   args$nclass, args$nclass_leaf,
                   is.init, init.param,
                   control$max.iter, control$tol,
                   control$verbose, control$per.iter)
   } else if (method == "vem") {
      fit <- emFit(data$y, args$nobs, args$nvar, args$ncat,
                   args$nlv, args$root - 1, args$leaf - 1,
                   args$u - 1, args$v - 1,
                   args$tree_index - 1, args$cstr_leaf - 1,
                   args$nclass, args$nclass_leaf,
                   is.init, init.param,
                   control$max.iter, control$tol,
                   control$verbose, control$per.iter)
   } else if (method == "gibbs") {
      fit <- emFit(data$y, args$nobs, args$nvar, args$ncat,
                   args$nlv, args$root - 1, args$leaf - 1,
                   args$u - 1, args$v - 1,
                   args$tree_index - 1, args$cstr_leaf - 1,
                   args$nclass, args$nclass_leaf,
                   is.init, init.param,
                   control$max.iter, control$tol,
                   control$verbose, control$per.iter)
   }

   logit <- logit_params(pi, tau, rho)

   loglik <- llFit(param, data$y, args$nobs, args$nvar, args$ncat,
                   args$nlv, args$root - 1, args$leaf - 1,
                   args$u - 1, args$v - 1,
                   args$tree_index - 1, args$cstr_leaf - 1,
                   args$nclass, args$nclass_leaf,)

   fit$param <- output_param(object$struct, args, fit$params)
   # fit$logit_param <- output_logit(object$struct, args, fit$params)
   # fit$se_param <- output_se(object$struct, args, fit$params)
   fit$posterior <- output_posterior(object$struct, args, data, fit$posterior)
   object$fit <- fit

   class(object) <- "catlvm.fit"
   object
}
