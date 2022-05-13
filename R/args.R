args_return <- function(struct) {
   constr <- struct$constr
   nvar <- sapply(struct$vars$manifest, length)

   args <- list(
      nlv = length(struct$label),
      root = as.numeric(struct$root),
      leaf = as.numeric(struct$leaf),
      u = as.numeric(struct$edge$child),
      v = as.numeric(struct$edge$parent),
      tree_index = as.numeric(struct$tree),
      cstr_leaf = constr$leaf,
      nclass = struct$nclass,
      nclass_leaf = constr$nclass_leaf,
      nroot = length(struct$root),
      nleaf = length(struct$leaf),
      nedge = nrow(struct$edge),
      nleaf_unique = length(constr$nclass_leaf),
      nvar = unname(nvar[!duplicated(constr$leaf)])
   )

   args
}

update_args <- function(args, data) {

   args$nobs <- data$nobs
   args$ncat <- unname(data$ncat[!duplicated(args$cstr_leaf)])

   npar_pi <- sum(args$nclass[args$root] - 1)
   npar_tau <- (args$nclass[args$u] - 1) %*% args$nclass[args$v]
   npar_rho <- sum(sapply(seq(args$nleaf_unique), function(v)
      (sum(args$ncat[[v]]) - args$nvar[v]) * args$nclass_leaf[v]))
   args$npar <- npar_pi + npar_tau + npar_rho

   args
}
