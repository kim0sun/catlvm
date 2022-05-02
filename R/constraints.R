proc_constr <- function(constraints, struct) {
   leaf <- struct$leaf
   edges <- struct$edges
   nclass <- struct$nclass
   leaf_constr <- letters[seq_len(length(leaf))]
   edge_constr <- letters[seq_len(nrow(edges))]

   for (i in seq_along(constraints)) {
      constr <- constraints[[i]]
      constr <- strsplit(gsub(" ", "", constr), "(->|~)")
      if (!all(unlist(constr) %in% struct$label)) next

      clen <- sapply(constr, length)
      if (all(clen == 1)) {
         leaf_constr[match(constr, leaf)] <- i
      }
      if (all(clen == 2)) {
         ind <- lapply(asplit(edges, 1), unname)
         edge_constr[match(constr, ind)] <- i
      }
   }

   nclass_leaf <- unlist(sapply(split(nclass[leaf], leaf_constr), unique))
   nclass_u <- sapply(split(nclass[edges$child], edge_constr), unique)
   if (is.list(nclass_u)) nclass_u <- numeric()
   nclass_v <- sapply(split(nclass[edges$parent], edge_constr), unique)
   if (is.list(nclass_v)) nclass_v <- numeric()

   list(leaf = as.numeric(factor(leaf_constr)),
        edge = as.numeric(factor(edge_constr)),
        nclass_leaf = unname(nclass_leaf),
        nclass_u = unname(nclass_u),
        nclass_v = unname(nclass_v))
}
