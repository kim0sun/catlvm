##' @export
print.catlvm <- function(x, ...) {
   cat("CATegorical Latent Variable Model\n")

   cat("\nLatent Variables (Root*) :")
   label <- x$model$label
   label[x$args$root] <- paste0(label[x$args$root], "*")
   mat <- rbind(label, x$model$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement Model:")
   vars <- x$model$vars$manifest
   formula <- sapply(vars, paste, collapse = ", ")
   ind <- letters[x$model$constr$cstr_leaf]
   mat = cbind(paste(names(vars), "-> {", formula, "}"), ind)
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)

   if (length(x$model$vars$latent) > 0) {
      cat("\nLatent Dependent structure:")
      ind <- LETTERS[x$model$constr$cstr_link]
      inds <- split(ind, x$model$links$parent, drop = TRUE)
      vars <- x$model$vars$latent
      formula <- sapply(vars, paste, collapse = ", ")
      constr <- sapply(inds, paste, collapse = ", ")
      mat = cbind(paste(names(vars), "-> {", formula, "}"), constr)
      dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
      print(mat[, -1, drop = FALSE], quote = FALSE)
   }
}

##' @export
print.catlvm.fit <- function(x, digits = 5, ...) {
   pi <- x$estimates$par$pi
   tau <- x$estimates$par$tau
   rho <- x$estimates$par$rho

   cat("CATegorical Latent Variable Model (estimated)\n")

   cat("\nRoot class prevalence (Pi) :\n")
   nr <- x$args$nroot
   nc <- x$model$nclass[x$args$root]
   mat <- matrix(NA, nr, max(nc))
   dimnames(mat) <- list(root = x$model$root, class = seq(max(nc)))
   for (r in seq_len(nr)) {
      mat[r, seq_len(nc[r])] <- pi[[r]]
   }
   print.table(round(mat, digits), quote = FALSE, ...)

   cat("\nMeasurement Model:")
   vars <- x$model$vars$manifest
   formula <- sapply(vars, paste, collapse = ", ")
   mat = cbind(paste(names(vars), "-> {", formula, "}"),
               letters[x$model$constr$cstr_leaf])
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)

   nlf <- x$args$nleaf_unique
   for (v in seq_len(nlf)) {
      cat("\n")
      cat(names(rho)[v])
      print.table(round(rho[[v]], digits), quote = FALSE, ...)
   }

   lstr <- x$model$vars$latent
   if (length(lstr) > 0) {
      cat("\nTransition path of latent variables :\n")

      for (d in seq_len(length(lstr))) {
         parent <- names(lstr)[d]
         nclass <- x$args$nclass[x$model$label == parent]
         lind <- which(x$model$links$parent == parent)
         child <- as.character(x$model$links$child[lind])
         path <- sapply(tau[x$args$cstr_link[lind]], apply, 2, which.max)
         mat <- cbind(seq(nclass), " ->", path)
         dimnames(mat) <- list(rep("", nrow(mat)), c(parent, "", child))
         print(mat, quote = FALSE, right = TRUE)
      }
   }
}


##' @export
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


# ##' @export
# summary.catlvm = function(object, ...) {
#
# }

# ##' @export
# coef.catlvm = function(object, ...) {
#
# }
# ##' @export
# vcov.catlvm = function(object, ...) {
#
# }
#
# ##' @export
# score.catlvm = function(x, ...) {
#
# }
# ##' @export
# reorder.catlvm = function(x, ...) {
#
# }
#
# ##' @export
# anova.catlvm = function(x, ...) {
#
# }
#
# ##' @export
# # predict
#
# ##' @export
# # posterior
#
# ##' @export
# # model.matrix
#
# ##' @export
# # model.frame
#
# ##' @export
# # confint
#
#
