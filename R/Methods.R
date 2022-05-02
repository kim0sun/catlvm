simulate.catlvm <- function(
   object, nsim = 500, ncat = NULL, pi = NULL, tau = NULL, rho = NULL
) {
   pi_valid <- function(pi, nclass) {
      if (is.matrix(pi)) {
         dimTF = dim(pi) == c(nclass)
         sumTF = all(colSums(pi) == 1)
         if (dimTF && sumTF) return(log(pi))
         else return(pi_gnr(nclass))
      } else if (is.numeric(pi)) {
         lenTF = length(pi) == c(nclass)
         sumTF = sum(pi) == 1
         if (lenTF && sumTF) return(matrix(log(pi), nclass))
         else return(pi_gnr(nclass))
      } else return(pi_gnr(nclass))
   }

   tau_valid <- function(tau, nk, nl) {
      if (is.matrix(tau)) {
         dimTF = dim(tau) == c(nk, nl)
         sumTF = all(colSums(tau) == 1)
         if (dimTF && sumTF) return(log(tau))
         else return(tau_gnr(nk, nl))
      } else return(tau_gnr(nk, nl))
   }

   rho_valid <- function(rho, nclass, ncat) {
      if (is.list(rho)) {
         dim1 = lapply(rho, dim)
         dim2 = lapply(ncat, function(m) c(m, nclass))
         dimTF = identical(dim1, dim2)
         sumTF = all(sapply(rho, function(x) all(colSums(x) == 1)))
         if (dimTF && sumTF) return(log(unlist(rho)))
         else return(rho_gnr(nclass, ncat))
      } else if (is.numeric(rho)) {
         ind = rep(rep(seq(length(ncat)), ncat), nclass)
         lrho = lapply(split(rho, ind), function(x)
            matrix(x, ncol = nclass))
         lenTF = length(rho) == nclass * sum(ncat)
         sumTF = all(sapply(lrho, function(x) all(colSums(x) == 1)))
         if (lenTF && sumTF) return(log(rho))
         else return(rho_gnr(nclass, ncat))
      } else return(rho_gnr(nclass, ncat))
   }

   nlv = object$args$nlv
   root = object$args$root
   leaf = object$args$leaf
   ulv = object$args$u
   vlv = object$args$v

   cstr_root = object$args$cstr_root
   cstr_leaf = object$args$cstr_leaf
   cstr_edge = object$args$cstr_edge

   nclass = object$args$nclass
   nclass_u = object$args$nclass_u
   nclass_v = object$args$nclass_v
   nclass_leaf = object$args$nclass_leaf

   nroot = object$args$nroot
   nleaf = object$args$nleaf
   nedge = object$args$nedge
   nleaf_unique = object$args$nleaf_unique
   nedge_unique = object$args$nedge_unique

   nvar = object$args$nvar
   ncat = if (is.null(ncat)) lapply(nvar, function(x) rep(2, x)) else ncat

   if (length(pi) < nroot)
      pi = lapply(seq_len(nroot), function(x) NULL)
   for (r in seq_len(nroot)) {
      pi[[r]] = pi_valid(pi[[r]], nclass[root[r]])
   }
   if (length(tau) < nedge_unique)
      tau = lapply(seq_len(nedge_unique), function(x) NULL)
   for (d in seq_len(nedge_unique)) {
      tau[[d]] = tau_valid(tau[[d]], nclass_u[d], nclass_v[d])
   }
   if (length(rho) < nleaf_unique)
      rho = lapply(seq_len(nleaf_unique), function(x) NULL)
   for (v in seq_len(nleaf_unique)) {
      rho[[v]] = rho_valid(rho[[v]], nclass_leaf[v], ncat[[v]])
   }

   ysim = ysim(nsim, nlv, root - 1, leaf - 1, ulv - 1, vlv - 1,
               nclass, nroot, nleaf, nedge, ncat,
               cstr_edge - 1, cstr_leaf - 1,
               pi, tau, rho, TRUE)

   args = list(nobs = nsim, nvar = nvar, ncat = ncat, nlv = nlv,
               root = root, leaf = leaf, ulv = ulv, vlv = vlv,
               cstr_root = cstr_root, cstr_leaf = cstr_leaf, cstr_edge = cstr_edge,
               nclass = nclass, nclass_leaf = nclass_leaf,
               nclass_u = nclass_u, nclass_v = nclass_v)

   list(response = ysim$y, class = ysim$class, args = args,
        params = list(pi = pi, tau = tau, rho = rho))
}


print.catlvm <- function(x, ...) {
   cat("CATegorical Latent Variable Model\n")

   cat("\nLatent Variables:")
   mat <- rbind(x$struct$label, x$struct$nclass)
   dimnames(mat) <- list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement Model:")
   vars <- x$struct$vars$manifest

   formula <- sapply(vars, paste, collapse = ", ")
   mat = cbind(paste(names(vars), "-> {", formula, "}"),
               letters[x$constraints$leaf])
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)

   edges <- apply(x$struct$edges, 1, paste, collapse = " -> ")
   if (length(x$struct$vars$latent) > 0) {
      cat("\nLatent Dependent Structure:")
      mat = cbind(edges, letters[x$constraints$edge])
      dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
      print(mat[, -1, drop = FALSE], quote = FALSE)
   }
}

estimate.catlvm = function(
   x, data = parent.frame(),
   method = "em", smoothing = TRUE,
   control = catlvm.control(), ...
) {

   class(res) = "catlvm"
}

# print
# logLik
# anova
# predict
# posterior
# reorder
# summary
# simulate
# coef
# vcov
# model.matrix
# model.frame
# confint

npar.catlvm = function(object, ...) {

}

nobs.catlvm = function(object, ...) {

}

summary.catlvm = function(object, ...) {

}

logLik.catlvm = function(object, ...) {

}

coef.catlvm = function(object, ...) {

}

vcov.catlvm = function(object, ...) {

}

score.catlvm = function(x, ...) {

}

catlvm.control = function(...) {
   list(maxiter = maxiter, tol = tol,
        init.param = init.param, nrep = nrep,
        na.action = na.action)
}
