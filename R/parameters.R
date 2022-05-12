##' @export
tau <- function(object, logit = FALSE) {
   if (logit) object$fit$logit_params$rho
   else object$fit$params$tau
}

##' @export
rho <- function(object, logit = FALSE) {
   if (logit) object$fit$logit_params$rho
   else object$fit$params$rho
}
