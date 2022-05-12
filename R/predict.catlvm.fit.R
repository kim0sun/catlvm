##' @export
predict.catlvm.fit <- function(object, variable = NULL, ...) {
   if (is.null(variable)) variable <- object$struct$label
   class <- sapply(variable, function(var) {
      post <- object$fit$posterior[[var]]
      apply(post, 1, which.max)
   })
   if (length(variable) == 1) return(class[,1])
   data.frame(class)
}



