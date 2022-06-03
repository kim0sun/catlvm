##' @export
print.catlvm <- function(x, digits = 5, ...) {
   cat("CATegorical Latent Variable Model\n")
   if (x$estimated) {
      pi <- x$estimates$par$pi
      tau <- x$estimates$par$tau
      rho <- x$estimates$par$rho
   }

   cat("\nLatent variables (Root*) :")
   label <- x$model$label
   label[x$args$root] <- paste0(label[x$args$root], "*")
   mat <- rbind(label, x$model$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   if (x$estimated) {
      cat("\nRoot class prevalence (Pi) :\n")
      nr <- x$args$nroot
      nc <- x$model$nclass[x$args$root]
      mat <- matrix(NA, nr, max(nc))
      dimnames(mat) <- list(root = x$model$root, class = seq(max(nc)))
      for (r in seq_len(nr)) {
         mat[r, seq_len(nc[r])] <- pi[[r]]
      }
      print.table(round(mat, digits), quote = FALSE, ...)
   }

   cat("\nMeasurement model:")
   vars <- x$model$vars$manifest
   formula <- sapply(vars, paste, collapse = ", ")
   ind <- letters[x$model$constr$cstr_leaf]
   mat = cbind(paste(names(vars), "-> {", formula, "}"), ind)
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)


   if (x$estimated) {
      nlf <- x$args$nleaf_unique
      for (v in seq_len(nlf)) {
         cat("\n")
         cat(names(rho)[v])
         print.table(round(rho[[v]], digits), quote = FALSE, ...)
      }
   }

   if (x$args$nlink > 0) {
      cat("\nLatent dependent structure:")
      vars <- x$model$vars$latent
      formula <- sapply(vars, paste, collapse = ", ")
      mat = cbind(paste(names(vars), "-> {", formula, "}"))
      dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
      print(mat[, -1, drop = FALSE], quote = FALSE)

      if (x$args$nlink > x$args$nlink_unique) {
         ind <- LETTERS[x$model$constr$cstr_link]
         parent <- split(x$model$links$parent, ind)
         child <- split(x$model$links$child, ind)
         mapping <- function(x, y) {
            ygrp <- split(y, x, drop = TRUE)
            yy <- lapply(ygrp, paste, collapse = ", ")
            paste0(names(yy), " -> ", yy, " ")
         }
         maps <- mapply(mapping, parent, child)
         mat <- matrix("", nrow = max(sapply(maps, length)), ncol = length(maps))
         for (i in seq_len(length(maps))) {
            mat[seq_len(length(maps[[i]])), i] = maps[[i]]
         }
         cat("\nDependency constraints:\n")
         dimnames(mat) = list(rep("", nrow(mat)), names(maps))
         print(mat, quote = FALSE)
         if (x$estimated) {
            cat("\nMost probable path:\n")
            path <- lapply(tau, apply, 2, which.max)
            mat <- sapply(path, function(x) paste(names(x), "->", x))
            dimnames(mat) = list(rep("", nrow(mat)), names(maps))
            print(mat, quote = FALSE)
         }
      }
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

   if (x$args$nlink > 0) {
      cat("\nLatent Dependent structure:")
      if (x$args$nlink == x$args$nlink_unique) {
         vars <- x$model$vars$latent
         formula <- sapply(vars, paste, collapse = ", ")
         mat <- cbind(paste(names(vars), "-> {", formula, "}"))
         dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
         print(mat[, -1, drop = FALSE], quote = FALSE)
         cat("\n")
         path <- lapply(tau, apply, 2, which.max)
         mat <- c()
         for (i in seq_len(length(tau))) {
            parent <- as.character(x$model$links[i,1])
            child <- as.character(x$model$links[i,2])
            mati <- matrix(" ", nrow = max(sapply(path, length)), 3)
            mati[seq_len(length(path[[i]])), 1] <- names(path[[i]])
            mati[seq_len(length(path[[i]])), 2] <- "->"
            mati[seq_len(length(path[[i]])), 3] <- path[[i]]
            colnames(mati) <- c(parent, "", child)
            mat <- cbind(mat, mati, "")
         }
         rownames(mat) = rep("", nrow(mat))
         print(mat, quote = FALSE, right = TRUE)
      } else {
         cat("\n")
         ind <- LETTERS[x$model$constr$cstr_link]
         parent <- split(x$model$links$parent, ind)
         child <- split(x$model$links$child, ind)
         mapping <- function(x, y) {
            ygrp <- split(y, x, drop = TRUE)
            yy <- lapply(ygrp, paste, collapse = ", ")
            paste0(names(yy), " -> ", yy, " ")
         }
         maps <- mapply(mapping, parent, child)
         mat <- matrix("", nrow = max(sapply(maps, length)), ncol = length(maps))
         for (i in seq_len(length(maps))) {
            mat[seq_len(length(maps[[i]])), i] = maps[[i]]
         }
         path <- lapply(tau, apply, 2, which.max)
         mat2 <- sapply(path, function(x) paste(names(x), "->", x))
         mat <- rbind(mat, "", mat2)
         dimnames(mat) = list(rep("", nrow(mat)), names(maps))
         print(mat, quote = FALSE)
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
