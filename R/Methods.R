simulate.catlvm = function(
   object, n = 500, pi = NULL, tau = NULL, rho = NULL
) {
   root = object$root
   leaf = object$leaf
   ulv = object$ulv
   vlv = object$vlv
   cstr_lf = object$cstr_lf
   cstr_lk = object$cstr_lk

   nc = object$nc
   unc = object$unc
   vnc = object$vnc
   nclf = object$nclf

   nleaf = object$nleaf
   nuleaf = object$nuleaf
   nlink = object$nlink
   nulink = object$nulink

   nvar = object$nvar
   ncat = object$ncat

   if (is.null(pi)) {
      pi = list()
      for (r in seq_along(root)) {
         pi[[r]] = exp(pi_gnr(nc[root[r]], 1))
      }
   }
   if (is.null(tau)) {
      tau = list()
      for (d in seq(nulink)) {
         tau[[d]] = exp(matrix(tau_gnr(vnc[d], unc[d], 1), vnc[d]))
      }
   }
   if (is.null(rho)) {
      rho = list()
      for (v in seq(nuleaf)) {
         rho[[v]] = exp(rho_gnr(nclf[v], ncat[[v]]))
      }
   }

   cls = list()
   for (r in seq_along(root)) {
      cls[[root[r]]] = root_gnr(n, nc[root[r]], pi[[r]])
   }
   for (d in rev(seq(ulv))) {
      u = ulv[d]; v = vlv[d]
      tau_d = tau[[cstr_lk[d]]]
      cls[[u]] = cls_gnr(n, nc[u], nc[v], cls[[v]], tau_d);
   }
   y = list()
   for (v in seq(nleaf)) {
      rho_v = rho[[cstr_lf[v]]]
      y[[v]] = y_gnr(n, nclf[cstr_lf[v]], ncat[[cstr_lf[v]]], cls[[leaf[v]]], rho_v)
   }

   y
}

catlvm = function(
   x = NULL, ..., data, subset, weight,
   constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula = x

   ltv <- list(label = c(), nclass = c())
   msr <- list(formula = c(), label = c(), manifest = list())
   str <- list(formula = c(), parent = c(), child = list())
   reg <- list(formula = c(), label = c(), covariates = list())
   cst <- list(leaf = list(), edge = list())

   rhs <- function(f) f[[length(f)]]
   lhs <- function(f) if (length(f) == 3) f[[2]] else  NULL
   nclass <- function(f) as.numeric(lhs(f)[[length(lhs(f))]])

   strf <- function(f) {
      vars = all.vars(f)
      formula = sapply(2:length(vars), function(i)
         paste(vars[1], "->", paste(vars[i])))
      list(formula = list(formula),
           parent = vars[1], child = list(vars[-1]))
   }
   msrf <- function(f) {
      vars = all.vars(f)
      list(formula = paste(vars[1], "->", paste(vars[-1], collapse = ", ")),
           label = vars[1], manifest = list(vars[-1]))
   }
   regf <- function(f) {
      list(formula = deparse(f))
   }
   typef <- function(f, lvs) {
      if (length(lhs(f)) > 2) {
         if (all(all.vars(rhs(f)) %in% lvs)) return("str")
         if (all(!(all.vars(rhs(f)) %in% lvs))) return("msr")
         else stop("Both manifest and latent variable are included in same formula.")
      } else return("reg")
   }
   labels <- sapply(formula, function(x) all.vars(lhs(x)))
   ltv = list(label = unique(labels),
              nclass = sapply(formula[!duplicated(labels)], function(x) nclass(x)))

   for (f in formula) {
      type = typef(f, ltv$label)
      if (type == "msr")
         msr = mapply(c, msr, msrf(f), SIMPLIFY = FALSE)
      if (type == "str")
         str = mapply(c, str, strf(f), SIMPLIFY = FALSE)
      if (type == "reg")
         reg = mapply(c, reg, regf(f), SIMPLIFY = FALSE)
   }

   msr_cstr = letters[seq(msr$formula)]
   str_cstr = letters[seq(unlist(str$formula))]
   for (i in seq_along(constraints)) {
      cstr = constraints[[i]]
      if (all(grepl("->", cstr)))
         str_cstr[which(unlist(str$formula) %in% cstr)] = i
      else {
         msr_cstr[which(msr$label %in% cstr)] = i
      }
   }
   msr$constraints = as.numeric(factor(msr_cstr))
   str$constraints = as.numeric(factor(str_cstr))

   root = setdiff(str$parent, unlist(str$child))
   leaf = msr$label

   res = list(estimated = FALSE)
   res$latentVariable = c(N = length(ltv$label), ltv,
                          root = root, leaf = leaf)
   res$measureModel = c(N = length(msr$formula), msr)
   res$latentStruct = c(N = length(str$formula), str)
   res$regressModel = c(N = length(reg$formula), reg)
   res$index = list(
      root <- sapply(root, function(x) which(ltv$label %in% x)),
      leaf <- sapply(leaf, function(x) which(ltv$label %in% x)),
      u <- sapply(unlist(str$child), function(x) which(ltv$label %in% x)),
      v <- sapply(rep(str$parent, sapply(str$child, length)),
                    function(x) which(ltv$label %in% x)),
      nclass_leaf = ltv$nclass[leaf[which(!duplicated(msr$constraints))]],
      nclass_u = ltv$nclass[ulv[which(!duplicated(str$constraints))]],
      nclass_v = ltv$nclass[vlv[which(!duplicated(str$constraints))]]
   )

   class(res) = "catlvm"
   res
}

f = list(l1[3] ~ x1 + x2,
         l2[3] ~ y1 + y2,
         l3[3] ~  z1 + z2,
         p[2] ~ l1 + l2 + l3,
         l1 ~ xx, p ~ yy)

a = catlvm(f, constraints = list(c("l1", "l2"), c("p -> l1", "p -> l2")))

print.catlvm <- function(x, ...) {
   cat("CATegorical Latent Variable Model",
       if (x$estimated) "(estimated)" ,"\n")

   cat("\nLatent Variables:")
   mat = rbind(x$latentVariable$label, x$latentVariable$nclass)
   dimnames(mat) = list(c(" Label:", "nclass:"), rep("", ncol(mat)))
   print(mat, quote = FALSE)

   cat("\nMeasurement Model:")
   msr = x$measureModel
   mat = cbind(msr$formula, letters[msr$constraints])
   dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
   print(mat[, -1, drop = FALSE], quote = FALSE)

   if (x$latentStruct$N > 0) {
      cat("\nLatent Dependent Structure:")
      str = x$latentStruct
      mat = cbind(unlist(str$formula), letters[str$constraints])
      dimnames(mat) = list(paste(" ", mat[,1], " "), rep("", ncol(mat)))
      print(mat[, -1, drop = FALSE], quote = FALSE)
   }

   if (x$regressModel$N > 0) {
      cat("\nRegression Model:")
      reg = x$regressModel
      mat = matrix(unlist(reg$formula), nrow = reg$N)
      rownames(mat) = paste(" ", unlist(reg$formula))
      print(mat[, NULL], quote = FALSE)
   }
}

estimate.catlvm = function(
   x, data = parent.frame(), subset, weights,
   method = "em", control = catlvm.control(), ...
) {

   class(res) = "catlvm"
}

npar.catlvm = function(object, ...) {

}

nobs.catlvm = function(object, ...) {

}

print.catlvm = function(x, ...) {

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
