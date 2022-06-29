#' @export
catlvm = function(
   x, ..., constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   model <- proc_formula(formula, constraints)

   res = list()
   res$model <- model
   res$args <- args_return(model)
   res$estimated <- FALSE

   class(res) <- "catlvm"
   res
}
