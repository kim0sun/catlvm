##' @export
print.catlvm <- function(x, ...) {
   cat("CATegorical Latent Variable Model\n")

   cat("\nLatent Variables (Root*) :")
   label <- x$struct$label
   label[x$args$root] <- paste0(label[x$args$root], "*")
   mat <- rbind(label, x$struct$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement Model:")
   vars <- x$struct$vars$manifest
   formula <- sapply(vars, paste, collapse = ", ")
   mat = cbind(paste(names(vars), "-> {", formula, "}"),
               letters[x$struct$constr$leaf_constr])
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)

   if (length(x$struct$vars$latent) > 0) {
      cat("\nLatent Dependent Structure:")
      vars <- x$struct$vars$latent
      formula <- sapply(vars, paste, collapse = ", ")
      mat = cbind(paste(names(vars), "-> {", formula, "}"))
      dimnames(mat) = list(paste(" ", mat[,1], " "), "")
      print(mat[, -1, drop = FALSE], quote = FALSE)
   }
}

##' @export
print.catlvm.fit <- function(x, ...) {
   pi <- x$fit$estimates$par$pi
   tau <- x$fit$estimates$par$tau
   rho <- x$fit$estimates$par$rho

   cat("\nCATegorical Latent Variable Model (estimated)\n")

   cat("\nRoot class prevalence (Pi) :\n")
   nr <- x$args$nroot
   nc <- x$struct$nclass[x$args$root]
   mat <- matrix(NA, nr, max(nc))
   dimnames(mat) <- list(root = x$struct$root, class = seq(max(nc)))
   for (r in seq_len(nr)) {
      mat[r, seq_len(nc[r])] <- pi[[r]]
   }
   print(mat, quote = FALSE, ...)

   cat("\nConditional transition probabilities (Tau) :")


   cat("\nItem response probabilities (Rho) :")

}


##' @export
summary.catlvm = function(object, ...) {

}

#' @export
print.catlvm.parameter <- function(x, ...) {
   y <- x
   attributes(y) <- attributes(x)[c("dim", "dimnames")]
   print.default(y)
}

##' @export
logLik.catlvm <- function(object, ...) {
   res <- if (is.null(obejct$fit)) NA
   else structure(
      object$fit$loglik,
      df = object$args$npar,
      nobs = object$data$nobs
   )
   class(res) <- "logLik"
   res
}


##' @export
coef.catlvm = function(object, ...) {

}
##' @export
vcov.catlvm = function(object, ...) {

}

##' @export
score.catlvm = function(x, ...) {

}
##' @export
reorder.catlvm = function(x, ...) {

}

##' @export
anova.catlvm = function(x, ...) {

}

##' @export
# predict

##' @export
# posterior

##' @export
# model.matrix

##' @export
# model.frame

##' @export
# confint







