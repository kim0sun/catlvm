##' @export
predict.catlvm.fit <- function(object, newdata = NULL, ...) {
   if (missing(newdata))
   calcll(unlist(object$estimates$logit))
}



