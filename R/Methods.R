simulate.catlvm = function(
   object, n = 500, ncat = NULL, pi = NULL, tau = NULL, rho = NULL
) {
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
   nclass_root = object$index$nclass_root
   nclass_leaf = object$index$nclass_leaf

   nleaf = object$measureModel$N
   nleaf_unique = length(object$index$unique_msr)
   nedge_unique = length(object$index$unique_str)

   nvar = object$index$nvar
   ncat = NULL
   ncat = if (is.null(ncat)) lapply(nvar, function(x) rep(2, x)) else ncat

   params = list()
   if (is.null(pi)) {
      pi = list()
      for (r in seq_along(root)) {
         pi[[r]] = exp(pi_gnr(nclass[root[r]], 1))
      }
      params$pi = pi
   }
   if (is.null(tau) && nedge_unique > 0) {
      tau = list()
      for (d in seq_len(nedge_unique)) {
         tau[[d]] = exp(matrix(tau_gnr(nclass_v[d], nclass_u[d], 1), nclass_v[d]))
      }
      params$tau = tau
   }
   if (is.null(rho)) {
      rho = list()
      for (v in seq_len(nleaf_unique)) {
         rho[[v]] = exp(rho_gnr(nclass_leaf[v], ncat[[v]]))
      }
      params$rho = rho
   }

   cls = list()
   for (r in seq_along(root)) {
      cls[[root[r]]] = root_gnr(n, nclass[root[r]], pi[[r]])
   }

   for (d in rev(seq_along(ulv))) {
      u = ulv[d]; v = vlv[d]
      tau_d = tau[[cstr_edge[d]]]
      cls[[u]] = cls_gnr(n, nclass[u], nclass[v], cls[[v]], tau_d);
   }
   y = list()
   for (v in seq_len(nleaf)) {
      rho_v = rho[[cstr_leaf[v]]]
      y[[v]] = y_gnr(n, nclass_leaf[cstr_leaf[v]], ncat[[cstr_leaf[v]]], cls[[leaf[v]]], rho_v)
   }

   list(y = y, nobs = n, nvar = nvar, ncat = ncat,
        nlv = length(cls), root = root, leaf = leaf,
        ulv = ulv, vlv = vlv,
        cstr_root = cstr_root, cstr_leaf = cstr_leaf, cstr_edge = cstr_edge,
        nclass = nclass, nclass_root = nclass_root, nclass_leaf = nclass_leaf,
        nclass_u = nclass_u, nclass_v = nclass_v,
        params = params)
}

catlvm = function(
   x = NULL, ..., data, subset, weight,
   constraints = NULL
) {
   if (!is.list(x)) formula <- list(x, ...)
   else formula = x

   rhs <- function(f) f[[length(f)]]
   lhs <- function(f) if (length(f) == 3) f[[2]] else  NULL
   ncls <- function(f) as.numeric(lhs(f)[[length(lhs(f))]])
   typef <- function(f, lvs) {
      if (length(lhs(f)) > 2) {
         if (all(all.vars(rhs(f)) %in% lvs)) return("str")
         if (all(!(all.vars(rhs(f)) %in% lvs))) return("msr")
         else stop("Both manifest and latent variable are included in same formula.")
      } else return("reg")
   }
   stage <- function(fs, lvs) {

      type = sapply(fs, function(f) typef(f, lvs))
      stages = integer(length(type))
      stages[which(type %in% c("msr", "reg"))] = 1


      all.vars(rhs(fs[[1]]))
   }

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
           label = vars[1], manifest = list(vars[-1]), nvar = length(vars) - 1)
   }
   regf <- function(f) {
      list(formula = deparse(f))
   }

   labels <- sapply(formula, function(x) all.vars(lhs(x)))
   terms  <- sapply(formula, function(x) all.vars(rhs(x)))
   if (any(duplicated(unlist(terms))))
      stop("More than 1 latent variable is measured by identical latent variable.")
   types  <- sapply(formula, typef, labels)
   stages <- integer(length(formula[types != "reg"]))
   names(stages) <- labels[types != "reg"]
   stages[which(types == "msr")] <- 1
   while (any(stages == 0)) {
      for (i in which(stages == 0)) {
         stage_terms <- stages[terms[[i]]]
         if (all(stage_terms > 0))
            stages[i] = max(stage_terms) + 1
      }
   }

   formula = c(formula[types != "reg"][order(stages)],
               formula[types == "reg"])
   labels <- sapply(formula, function(x) all.vars(lhs(x)))
   nclass <- sapply(formula[!duplicated(labels)], ncls)
   ltv <- list(label = unique(labels), nclass = nclass)
   msr <- list(formula = c(), label = c(), manifest = list(), nvar = c())
   str <- list(formula = c(), parent = c(), child = list())
   reg <- list(formula = c(), label = c(), covariates = list())

   for (f in formula) {
      type = typef(f, ltv$label)
      if (type == "msr")
         msr = mapply(c, msr, msrf(f), SIMPLIFY = FALSE)
      if (type == "str")
         str = mapply(c, str, strf(f), SIMPLIFY = FALSE)
      if (type == "reg")
         reg = mapply(c, reg, regf(f), SIMPLIFY = FALSE)
   }

   leaf = msr$label
   root = if (length(str$formula) > 0) {
      setdiff(str$parent, unlist(str$child))
   } else msr$label

   ltv_cstr = letters[seq_along(root)]
   msr_cstr = letters[seq_along(msr$formula)]
   str_cstr = letters[seq_along(unlist(str$formula))]
   if (is.null(constraints)) {
      msr$constraints = seq_len(length(msr$label))
      str$constraints = seq_len(length(str$label))
      ltv$constraints = seq_len(length(root))
   } else {
      for (i in seq_along(constraints)) {
         cstr = constraints[[i]]
         if (all(grepl("->", cstr)))
            str_cstr[which(unlist(str$formula) %in% cstr)] = i
         else if (all(cstr %in% msr$label)) {
            msr_cstr[which(msr$label %in% cstr)] = i
         } else {
            ltv_cstr[which(root %in% cstr)] = i
         }
      }
      msr$constraints = as.numeric(factor(msr_cstr))
      str$constraints = as.numeric(factor(str_cstr))
      ltv$constraints = as.numeric(factor(ltv_cstr))
   }

   res = list(estimated = FALSE)
   res$latentVariable = c(N = length(ltv$label), ltv,
                          root = root, leaf = leaf)
   res$measureModel = c(N = length(msr$formula), msr)
   res$latentStruct = c(N = length(str$formula), str)
   res$regressModel = c(N = length(reg$formula), reg)

   root = sapply(root, function(x) which(ltv$label %in% x))
   leaf = sapply(leaf, function(x) which(ltv$label %in% x))
   u = unlist(sapply(unlist(str$child), function(x) which(ltv$label %in% x)))
   v = unlist(sapply(rep(str$parent, sapply(str$child, length)),
              function(x) which(ltv$label %in% x)))
   unique_rt  = which(!duplicated(ltv$constraints))
   unique_msr = which(!duplicated(msr$constraints))
   unique_str = which(!duplicated(str$constraints))

   res$index = list(
      root = root, leaf = leaf, u = u, v = v,
      nvar = msr$nvar[unique_msr],
      unique_rt = unique_rt,
      unique_msr = unique_msr,
      unique_str = unique_str,
      nclass_root = ltv$nclass[root[unique_rt]],
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


f = list(L1[3] ~ x1 + x2,
         L2[3] ~ y1 + y2,
         L3[3] ~  z1 + z2,
         M1[3] ~ x12 + x22,
         M2[3] ~ y12 + y22,
         M3[3] ~  z12 + z22,
         N1[3] ~ x13 + x23,
         N2[3] ~ y13 + y23,
         N3[3] ~  z13 + z23,
         P1[3] ~ L1 + L2 + L3,
         P2[3] ~ M1 + M2 + M3,
         P3[3] ~ N1 + N2 + N3,
         U[2] ~ P1 + P2 + P3)
lvm = catlvm(f, constraints = list(c("L1", "L2", "L3"), c("M1", "M2", "M3"), c("N1", "N2", "N3")))
lvm = catlvm(L1[3] ~ X1 + X2 + X3, L2[3] ~ Y1 + Y2 + Y3)
lvm
y = simulate(lvm, 50)
fit = treeFit(unlist(y$y), y$nobs, y$nvar, y$ncat,
        y$nlv, y$root - 1, y$leaf - 1, y$ulv - 1, y$vlv - 1,
        y$cstr_leaf - 1, y$cstr_edge - 1,
        y$nclass, y$nclass_leaf,
        y$nclass_u, y$nclass_v, 10, 1e-3)



exp(fit$logliks)
log(c(exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(2, 8, 13)])),
   exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(4, 10, 15)])),
   exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(6, 12, 17)]))))
log(exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(2, 7, 14)])) +
      exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(4, 9, 16)])) +
      exp(fit$pi[[1]][1] + sum((fit$rho[[1]])[c(6, 11, 18)])))
log(sum(exp(fit$a[[1]] + fit$b[[1]])[1:3]))
(fit$logliks)
fit$pi[[1]] +
fit$logliks

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
