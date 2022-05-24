##' @export
catlvm = function(
   x, ..., data, constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   struct <- proc_formula(formula, constraints)

   res = list()
   res$struct <- struct
   res$args <- args_return(struct)
   if (!missing(data)) {
      data <- proc_data(data, struct)
      res$data <- data
      res$args <- update_args(res$args, data)
   }

   class(res) <- "catlvm"
   res
}
