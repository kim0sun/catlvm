##' @export
catlvm = function(
   x, ..., data, group, constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   model <- proc_formula(formula, constraints)

   res = list()
   res$model <- model
   res$args <- args_return(model)
   if (!missing(data)) {
      data <- proc_data(data, model)
      res$data <- data
      res$args <- update_args(res$args, data)
   }
   res$estimated <- FALSE

   class(res) <- "catlvm"
   res
}
