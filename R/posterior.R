##' @export
posterior <- function(object, ...) UseMethod("posterior")

##' @export
posterior.catlvm.fit <- function(object, variable = NULL, digits = 3, ...) {
   post <- object$posterior

   if (is.null(variable)) variable <- object$struct$root[1]

   if (length(variable) == 1) out <- post[[variable]]
   else out <- post[variable]

   class(out) <- "posterior.catlvm"
   out
}

##' @export
print.posterior.catlvm <- function(x, n = 6L, ...) {
   if (is.list(x)) print(lapply(x, head, n, ...))
   else print(head(x, n, ...))
   invisible()
}


