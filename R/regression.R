#'
#' @export
regress <- function(object, ...) UseMethod("regress")

#'
#' @export
regress.catlvm <- function(
   object, formula, data = parent.frame(),
   imputation = c("modal", "prob"),
   method = c("3naive", "3BCH", "3ML"), ...
) {
   if (!object$estimated) stop("Latent variable model should be estimated.")

   # Import
   labels <- all.vars(formula)
   latent <- labels[labels %in% object$model$label]

   imputed <- sapply(object$posterior[latent], impute, imputation)
   data <- cbind(data, imputed)
   mf <- model.frame(formula, data)

   # naive (biased)
   glm(formula, mf, family = "binomial")

   # BCH


   # ML

   return(object)
}

impute <- function(posterior, imputation) switch(
   imputation, modal = apply(posterior, 1, which.max),
   prob = apply(posterior, 1, function(x)
      sample(seq_len(length(x)), 1, prob = x))
)
