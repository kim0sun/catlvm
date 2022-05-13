##' @export
catlvm = function(
   x = NULL, ..., data = NULL, constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   struct <- proc_formula(formula, constraints)

   res = list()
   res$struct <- struct
   res$args <- args_return(struct)
   if (!is.null(data)) {
      data <- proc_data(data, struct)
      rest$data <- data
      res$args <- update_args(args, data)
   }

   class(res) <- "catlvm"
   res
}
