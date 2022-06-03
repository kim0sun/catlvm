##' @export
posterior <- function(object, ...) UseMethod("posterior")

##' @export
posterior.catlvm.fit <- function(object, variable = NULL, ...) {
   post <- object$posterior

   if (is.null(variable)) variable <- object$model$root[1]

   if (length(variable) == 1) {
      out <- post[[variable]]
      names(dimnames(out))[2] <- paste0(variable, ": class")
   }
   else out <- post[variable]

   class(out) <- "posterior"
   out
}

##' @export
print.posterior <- function(x, n = 6L, ...) {
   if (is.list(x)) print(lapply(x, head, n, ...))
   else print(head(x, n, ...))
   invisible()
}
