##' @export
posterior <- function(object, ...) UseMethod("posterior")

##' @export
posterior.catlvm.fit <- function(object, variable = NULL, digits = 3, ...) {
   post <- object$fit$posterior

   if (is.null(variable)) variable <- object$struct$root[1]

   post[[variable]]
}

