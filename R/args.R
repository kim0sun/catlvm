args_return <- function(model) {
   constr <- model$constr
   nvar <- sapply(model$vars$manifest, length)

   args <- list(
      nlv = length(model$label),
      root = as.numeric(model$root),
      leaf = as.numeric(model$leaf),
      u = as.numeric(model$links$child),
      v = as.numeric(model$links$parent),
      tree_index = as.numeric(model$tree),
      cstr_leaf = constr$cstr_leaf,
      cstr_link = constr$cstr_link,
      nclass = model$nclass,
      nclass_root = model$nclass[as.numeric(model$root)],
      nclass_leaf = constr$nclass_leaf,
      nclass_u = constr$nclass_u,
      nclass_v = constr$nclass_v,
      nroot = length(model$root),
      nleaf = length(model$leaf),
      nlink = nrow(model$links),
      nleaf_unique = length(constr$nclass_leaf),
      nlink_unique = length(unique(constr$cstr_link)),
      nvar = unname(nvar[!duplicated(constr$cstr_leaf)])
   )

   args
}

update_args <- function(args, data) {
   args$nobs <- data$nobs
   ncat <- data$ncat[!duplicated(args$cstr_leaf)]
   names(ncat) <- letters[unique(args$cstr_leaf)]
   args$ncat <- ncat

   npar_pi <- sum(args$nclass[args$root] - 1)
   npar_tau <- c((args$nclass_u - 1) %*% args$nclass_v)
   npar_rho <- sum(sapply(seq(args$nleaf_unique), function(v)
      (sum(args$ncat[[v]]) - args$nvar[v]) * args$nclass_leaf[v]))
   args$npar <- c(pi = npar_pi, tau = npar_tau, rho = npar_rho)

   args
}
