args_return <- function(model) {
   constr <- model$constr
   nvar <- sapply(model$vars$manifest, length)

   args <- list(
      nlv = length(model$label),
      root = as.numeric(model$root),
      leaf = as.numeric(model$leaf),
      u = as.numeric(model$edge$child),
      v = as.numeric(model$edge$parent),
      tree_index = as.numeric(model$tree),
      cstr_leaf = constr$leaf,
      nclass = model$nclass,
      nclass_leaf = constr$nclass_leaf,
      nroot = length(model$root),
      nleaf = length(model$leaf),
      nedge = nrow(model$edge),
      nleaf_unique = length(constr$nclass_leaf),
      nvar = unname(nvar[!duplicated(constr$leaf)])
   )

   args
}

update_args <- function(args, data) {
   args$nobs <- data$nobs
   ncat <- data$ncat[!duplicated(args$cstr_leaf)]
   names(ncat) <- letters[unique(args$cstr_leaf)]
   args$ncat <- ncat

   npar_pi <- sum(args$nclass[args$root] - 1)
   npar_tau <- c((args$nclass[args$u] - 1) %*% args$nclass[args$v])
   npar_rho <- sum(sapply(seq(args$nleaf_unique), function(v)
      (sum(args$ncat[[v]]) - args$nvar[v]) * args$nclass_leaf[v]))
   args$npar <- c(pi = npar_pi, tau = npar_tau, rho = npar_rho)

   args
}
