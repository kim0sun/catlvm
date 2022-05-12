##' @export
catlvm = function(
   x = NULL, ..., data = NULL, constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   struct <- proc_formula(formula, constraints)
   data <- proc_data(data, struct)

   res = list()
   res$struct <- struct
   res$data <- data
   res$args <- args_return(struct, data)

   class(res) <- "catlvm"
   res
}
