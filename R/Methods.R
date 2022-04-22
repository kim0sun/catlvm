simulate.catlvm <- function(
   object, nsim = 500, ncat = NULL, pi = NULL, tau = NULL, rho = NULL
) {
   pi_valid <- function(pi, nobs, nclass) {
      if (is.matrix(pi)) {
         dimTF = dim(pi) == c(nclass, nobs)
         sumTF = all(colSums(pi) == 1)
         if (dimTF && sumTF) return(log(pi))
         else return(pi_gnr(nclass, nobs))
      } else if (is.numeric(pi)) {
         lenTF = length(pi) == c(nclass)
         sumTF = sum(pi) == 1
         if (lenTF && sumTF) return(matrix(log(pi), nclass, nobs))
         else return(pi_gnr(nclass, nobs))
      } else return(pi_gnr(nclass, nobs))
   }

   tau_valid <- function(tau, nobs, nk, nl) {
      if (is.matrix(tau)) {
         dimTF1 = dim(pi) == c(nk, nl)
         dimTF2 = dim(pi) == c(nk * nl, nobs)
         sumTF = all(colSums(pi) == 1)
         if (dimTF1 && sumTF) return(matrix(log(tau), nk * nl, nobs))
         else if (dimTF2 && sumTF) return(log(tau))
         else return(tau_gnr(nk, nl, nobs))
      } else return(tau_gnr(nk, nl, nobs))
   }

   rho_valid <- function(rho, nobs, nclass, ncat) {
      if (is.list(rho)) {
         dim1 = lapply(rho, dim)
         dim2 = lapply(ncat, function(m) c(m, nclass))
         dimTF = identical(dim1, dim2)
         sumTF = all(sapply(rho, function(x) all(colSums(x) == 1)))
         if (dimTF && sumTF) return(log(unlist(rho)))
         else return(rho_gnr(nclass, ncat))
      } else if (is.numeric(rho)) {
         ind = rep(seq(length(ncat)), nclass * ncat)
         lrho = lapply(split(rho, ind), function(x)
            matrix(x, nrow = nclass))
         lenTF = length(rho) == nclass * sum(ncat)
         sumTF = all(sapply(lrho, function(x) all(colSums(x) == 1)))
         if (lenTF && sumTF) return(log(rho))
         else return(rho_gnr(nclass, ncat))
      } else return(rho_gnr(nclass, ncat))
   }

   nlv = object$latentVariable$N
   root = object$index$root
   leaf = object$index$leaf
   ulv = object$index$u
   vlv = object$index$v
   cstr_root = object$latentVariable$constraints
   cstr_leaf = object$measureModel$constraints
   cstr_edge = object$latentStruct$constraints

   nclass = object$latentVariable$nclass
   nclass_u = object$index$nclass_u
   nclass_v = object$index$nclass_v
   nclass_leaf = object$index$nclass_leaf

   nroot = length(root)
   nleaf = object$measureModel$N
   nleaf_unique = length(object$index$unique_msr)
   nedge_unique = length(object$index$unique_str)

   nvar = object$index$nvar
   ncat = if (is.null(ncat)) lapply(nvar, function(x) rep(2, x)) else ncat

   params = list()
   params$pi = list()
   if (length(pi) < nroot)
      pi = lapply(seq_len(nroot), function(x) NULL)
   for (r in seq_along(root)) {
      params$pi[[r]] = pi_valid(pi[[r]], nsim, nclass[root[r]])
   }
   params$tau = list()
   if (length(tau) < nedge_unique)
      tau = lapply(seq_len(nedge_unique), function(x) NULL)
   for (d in seq_len(nedge_unique)) {
      params$tau[[d]] = tau_valid(tau[[d]], nsim, nclass_u[d], nclass_v[d])
   }
   params$rho = list()
   if (length(rho) < nleaf_unique)
      pi = lapply(seq_len(nleaf_unique), function(x) NULL)
   for (v in seq_len(nleaf_unique)) {
      params$rho[[v]] = rho_valid(rho[[v]], nobs, nclass_leaf[v], ncat[[v]])
   }

   cls = list()
   for (r in seq_along(root)) {
      cls[[root[r]]] = root_gnr(nsim, nclass[root[r]], params$pi[[r]])
   }

   for (d in rev(seq_along(ulv))) {
      u = ulv[d]; v = vlv[d]
      tau_d = params$tau[[cstr_edge[d]]]
      cls[[u]] = cls_gnr(nsim, nclass[u], nclass[v], cls[[v]], tau_d);
   }
   y = list()
   for (v in seq_len(nleaf)) {
      u = leaf[v]
      cstr = cstr_leaf[v]
      rho_v = params$rho[[cstr]]
      y[[v]] = y_gnr(nsim, nclass[u], ncat[[cstr]], cls[[u]], rho_v)
   }

   index = list(
      nobs = nsim, nvar = nvar, ncat = ncat, nlv = nlv,
      root = root, leaf = leaf, ulv = ulv, vlv = vlv,
      cstr_leaf = cstr_leaf, cstr_edge = cstr_edge,
      nclass = nclass, nclass_leaf = nclass_leaf,
      nclass_u = nclass_u, nclass_v = nclass_v
   )

   list(response = y, class = cls, index = index, params = params)
}

catlvm = function(
   x = NULL, ..., data, subset, weight,
   constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula <- x

   rhs <- function(f) f[[length(f)]]
   lhs <- function(f) if (length(f) == 3) f[[2]] else  NULL
   ncls <- function(f) as.numeric(lhs(f)[[length(lhs(f))]])
   typef <- function(f, lvs) {
      ll = length(lhs(f)) > 2
      rl = all(all.vars(rhs(f)) %in% lvs)
      rm = all(!(all.vars(rhs(f)) %in% lvs))
      if (rl) return(2)
      if (ll && rm) return(1)
      if (!ll && rm) return(3)
      if (!rl && !rm) {
         stop("Formula wrong.")
      }
   }

   strf <- function(f) {
      vars = all.vars(f)
      f = sapply(2:length(vars), function(i)
         paste(vars[1], "->", paste(vars[i])))
      list(formula = f, parent = rep(vars[1], length(vars[-1])),
           child = vars[-1])
   }
   msrf <- function(f) {
      vars = all.vars(f)
      f = paste(vars[1], "->", paste(vars[-1], collapse = ", "))
      lb = vars[1]
      m = list(vars[-1])
      names(m) = lb
      list(formula = f, label = lb,
           manifest = m, nvar = length(vars[-1]))
   }
   regf <- function(f) {
      list(formula = deparse(f))
   }

   labels <- sapply(formula, function(x) all.vars(lhs(x)))
   terms  <- lapply(formula, function(x) all.vars(rhs(x)))
   if (any(duplicated(unlist(terms))))
      stop("Some variables measures more than 1 variable.")
   types  <- sapply(formula, typef, labels)

   stages <- rep(-1, sum(types < 3))
   names(stages) <- labels[types < 3]
   stages[which(types == 1)] <- 0
   iter = 0
   while (any(stages < 0)) {
      iter = iter + 1
      unstaged = which(stages < 0)
      bot = sapply(terms[unstaged], function(x)
         all(!x %in% labels[unstaged]))
      stages[unstaged[bot]] = iter
   }

   formula = c(formula[types < 3][order(stages)],
               formula[types == 3])
   labels <- sapply(formula, function(x) all.vars(lhs(x)))
   nclass <- sapply(formula[!duplicated(labels)], ncls)

   msr <- list(formula = c(), label = c(), manifest = list(), nvar = c())
   str <- list(formula = c(), parent = c(), child = c())
   reg <- list(formula = c(), label = c(), covariates = list())

   for (f in formula) {
      type = typef(f, unique(labels))
      if (type == 1)
         msr = mapply(c, msr, msrf(f), SIMPLIFY = FALSE)
      if (type == 2)
         str = mapply(c, str, strf(f), SIMPLIFY = FALSE)
      if (type == 3)
         reg = mapply(c, reg, regf(f), SIMPLIFY = FALSE)
   }

   leaf_lab = msr$label
   root_lab = setdiff(sapply(formula, function(x) all.vars(lhs(x))),
                      unlist(sapply(formula, function(x) all.vars(rhs(x)))))
   ltv <- list(label = unique(labels), nclass = nclass,
               root = root_lab, leaf = leaf_lab)

   msr_cstr = letters[seq_along(msr$formula)]
   str_cstr = letters[seq_along(unlist(str$formula))]
   if (is.null(constraints)) {
      msr$constraints = seq_len(length(msr$label))
      str$constraints = seq_len(length(str$formula))
   } else {
      for (i in seq_along(constraints)) {
         cstr = constraints[[i]]
         if (all(grepl("->", cstr)))
            str_cstr[which(unlist(str$formula) %in% cstr)] = i
         else if (all(cstr %in% msr$label)) {
            msr_cstr[which(msr$label %in% cstr)] = i
         }
      }
      msr$constraints = as.numeric(factor(msr_cstr))
      str$constraints = as.numeric(factor(str_cstr))
   }

   res = list(estimated = FALSE)
   res$latentVariable = c(N = length(ltv$label), ltv)
   res$measureModel = c(N = length(msr$formula), msr)
   res$latentStruct = c(N = length(str$formula), str)
   res$regressModel = c(N = length(reg$formula), reg)

   root = sapply(root_lab, function(x) which(ltv$label %in% x))
   leaf = sapply(leaf_lab, function(x) which(ltv$label %in% x))
   u = unlist(sapply(unlist(str$child), function(x) which(ltv$label %in% x)))
   v = unlist(sapply(rep(str$parent, sapply(str$child, length)),
              function(x) which(ltv$label %in% x)))
   unique_rt  = which(!duplicated(ltv$constraints))
   unique_msr = which(!duplicated(msr$constraints))
   unique_str = which(!duplicated(str$constraints))

   res$index = list(
      root = root, leaf = leaf, u = u, v = v,
      nvar = msr$nvar[unique_msr],
      unique_msr = unique_msr,
      unique_str = unique_str,
      nclass_leaf = ltv$nclass[leaf[unique_msr]],
      nclass_u = ltv$nclass[u[unique_str]],
      nclass_v = ltv$nclass[v[unique_str]]
   )

   class(res) = "catlvm"
   res
}

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
lca   = catlvm(L1[2] ~ X1 + X2 + X3)
lcas  = catlvm(L1[2] ~ X1 + X2 + X3, L2[2] ~ Y1 + Y2 + Y3, L3[2] ~ Z1 + Z2 + Z3)
jlca  = catlvm(L1[2] ~ X1 + X2 + X3, L2[2] ~ Y1 + Y2 + Y3, L3[2] ~ Z1 + Z2 + Z3, L1 ~ L2, L2 ~ L3)
lcpa  = catlvm(L1[2] ~ X1 + X2 + X3, L2[2] ~ Y1 + Y2 + Y3, L3[2] ~ Z1 + Z2 + Z3, P1[2] ~ L1 + L2 + L3,
              constraints = list(c("L1", "L2", "L3")))
jlcpa = catlvm(L1[2] ~ X11 + X21 + X31, M1[2] ~ Y11 + Y21 + Y31, N1[2] ~ Z11 + Z21 + Z31,
               L2[2] ~ X12 + X22 + X32, M2[2] ~ Y12 + Y22 + Y32, N2[2] ~ Z12 + Z22 + Z32,
               L3[2] ~ X13 + X23 + X33, M3[2] ~ Y13 + Y23 + Y33, N3[2] ~ Z13 + Z23 + Z33,
               J1[3] ~ L1 + M1 + N1, J2[3] ~ L2 + M2 + N2, J3[3] ~ L3 + M3 + N3,
               JP[3] ~ J1 + J2 + J3,
               constraints = list(c("L1", "L2", "L3"), c("M1", "M2", "M3"), c("N1", "N2", "N3")))
lta = catlvm(L1[2] ~ X11 + X21 + X31, L2[2] ~ X12 + X22 + X32, L3[2] ~ X13 + X23 + X33, L1 ~ L2, L2 ~ L3,
             constraints = list(c("L1", "L2", "L3")))

lvm = lcas
y = simulate(lvm, 50)
fit2 = treeFit(unlist(y$response), y$index$nobs, y$index$nvar, y$index$ncat, y$index$nlv,
              y$index$root - 1, y$index$leaf - 1, y$index$ulv - 1, y$index$vlv - 1,
              y$index$cstr_leaf - 1, y$index$cstr_edge - 1,
              y$index$nclass, y$index$nclass_leaf, y$index$nclass_u, y$index$nclass_v,
              rep(TRUE, 3), y$params, 5000, 1e-5)

fit1$loglik
fit2$loglik
t(sapply(fit1$param$rho, function(x) round(exp(x), 3)))
t(sapply(fit2$param$rho, function(x) round(exp(x), 3)))

lapply(y$params$pi, function(x) round(exp(x), 3))
lapply(fit$pi, function(x) round(exp(x), 3))
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
