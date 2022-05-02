args_return <- function(res) {
   struct = res$struct
   constr = res$constraints
   # data_attr = res$data

   list(
      nlv = length(struct$label),
      root = as.numeric(struct$root),
      leaf = as.numeric(struct$leaf),
      u = as.numeric(struct$edge$child),
      v = as.numeric(struct$edge$parent),
      cstr_root = as.numeric(struct$tree),
      cstr_leaf = constr$leaf,
      cstr_edge = constr$edge,
      nclass = struct$nclass,
      nclass_u = constr$nclass_u,
      nclass_v = constr$nclass_v,
      nclass_leaf = constr$nclass_leaf,
      nroot = length(struct$root),
      nleaf = length(struct$leaf),
      nedge = nrow(struct$edge),
      nleaf_unique = length(constr$nclass_leaf),
      nedge_unique = length(constr$nclass_u),
      nvar = sapply(struct$vars$manifest, length)
   )
}
