catlvm = function(
   x = NULL, ..., data, subset, weight,
   constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   latent_struct <- proc_formula(formula)
   latent_constr <- proc_constr(constraints, struct)
   # data_attributes <- proc_data(data, latent_struct)

   res = list()
   res$struct <- latent_struct
   res$constraints <- latent_constr
   # res$data <- data_attributes
   res$args <- args_return(res)

   class(res) = "catlvm"
   res
}
