#' @export
predict.catlvm <- function(object, newdata, label, type = c("class", "posterior")) {
   type <- match.arg(type)
   if (!object$estimated) stop("Latent variable model should be estimated.")

   if (missing(label))
       label <- object$model$label

   if (missing(newdata))
      posterior = object$posterior
   else {
      data <- proc_data(newdata, object$model, FALSE)
      args <- object$args

      post <- calcPost(
         args$log_par, data$y, data$nobs, args$nvar, args$ncat,
         args$nlv, args$nroot, args$nlink, args$nleaf,
         args$nlink_unique, args$nleaf_unique,
         args$root - 1, args$tree_index - 1,
         args$u - 1, args$v - 1, args$leaf - 1,
         args$cstr_link - 1, args$cstr_leaf - 1,
         args$nclass, args$nclass_leaf,
         args$nclass_u, args$nclass_v
      )
      posterior <- output_posterior(post[label], object$model, data)
   }

   switch(type, class = if (length(label) == 1) impute(posterior[[label]])
      else sapply(posterior[label], impute),
      posterior = if (length(label) == 1) posterior[[label]]
      else posterior[label])
}
