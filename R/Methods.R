#' @export
print.catlvm <- function(x, digits = 5, ...) {
   cat("CATegorical Latent Variable Model\n")
   cat("\nLatent variables (Root*) :")
   label <- x$model$label
   label[x$args$root] <- paste0(label[x$args$root], "*")
   mat <- rbind(label, x$model$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement model:\n")
   vars <- x$model$vars$manifest
   collapsed <- sapply(vars, paste, collapse = ", ")
   formula <- paste(names(vars), "-> {", collapsed, "}")
   constr <- letters[x$model$constr$cstr_leaf]
   mat <- rbind(cbind(formula, " ", constr), "")
   dimnames(mat) <- list(rep("", nrow(mat)), mat[1,])
   print(mat[-1, , drop = FALSE], quote = FALSE)

   if (x$args$nlink > 0) {
      cat("Latent dependent structure:\n")
      vars <- x$model$vars$latent
      formula <- sapply(vars, paste, collapse = ", ")
      fmlae <- paste(names(vars), "-> {", formula, "}")
      cat("", paste(fmlae, collapse = "\n "), "\n")

      cstr <- LETTERS[x$model$constr$cstr_link]
      parent <- split(x$model$links$parent, cstr)
      child <- split(x$model$links$child, cstr)
      mapping <- function(x, y) paste(x, "->", y)
      maps <- mapply(mapping, parent, child, SIMPLIFY = FALSE)
      mat <- matrix("", nrow = max(sapply(maps, length)), ncol = length(maps))
      for (i in seq_len(length(maps))) {
         mat[seq_len(length(maps[[i]])), i] = maps[[i]]
      }
      cat("\nDependency constraints:\n")
      dimnames(mat) = list(rep("", nrow(mat)), names(maps))
      print(mat, quote = FALSE)
   }
}

#' @export
logLik.catlvm <- function(object, ...) {
   res <- if (is.null(obejct$estimate)) NA
   else structure(
      sum(object$loglik),
      df = object$args$npar,
      nobs = object$data$nobs
   )
   class(res) <- "logLik"
   res
}

# #' @export
# coef.catlvm = function(object, ...) {
#
# }
# #' @export
# vcov.catlvm = function(object, ...) {
#
# }
#
# #' @export
# score.catlvm = function(x, ...) {
#
# }
# #' @export
# reorder.catlvm = function(x, ...) {
#
# }
#
# #' @export
# anova.catlvm = function(x, ...) {
#
# }
#
# #' @export
# # predict
#
# #' @export
# # posterior
#
# #' @export
# # model.matrix
#
# #' @export
# # model.frame
#
# #' @export
# # confint
#
#
