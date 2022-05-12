label_formula <- function(fmla) {
   if (length(fmla) == 2) {
      return(NULL)
   }
   all.vars(fmla)[1]
}

nclass_formula <- function(fmla) {
   lhs <- fmla[[2]]
   part_lhs <- as.character(lhs)
   nclass <- setdiff(part_lhs, all.names(lhs))
   return(as.numeric(nclass))
}

get_lformula <- function(fmla) {
   if (!inherits(fmla, "formula"))
      fmla <- as.formula(fmla)

   attr(fmla, "label") = label_formula(fmla)
   attr(fmla, "nclass") = nclass_formula(fmla)
   attr(fmla, "vars") = labels(terms(fmla))

   fmla
}

proc_edges <- function(init_edges) {
   label <- levels(init_edges$ind)
   edges <- init_edges
   index_root <- logical(nrow(edges))
   stage <- integer(nrow(edges))
   while (any(!index_root)) {
      copy_edges <- edges
      names(copy_edges) <- c("ind", "ind2")
      merged <- merge(edges, copy_edges, all.x = TRUE)
      rooted <- which(is.na(merged[[3]]))
      stage[rooted] <- stage[rooted] + 1
      index_root[rooted] <- TRUE
      merged$ind2[rooted] <- merged$ind[rooted]
      edges <- data.frame(
         values = merged$values,
         ind = merged$ind2
      )
   }

   rt <- data.frame(root = edges$ind, lv = init_edges$ind)
   root <- unique(rt$root)
   leaf <- unique(init_edges$ind[!init_edges$values %in% label])
   tree <- droplevels(unique(rt[order(rt$lv),])$root)
   edges <- init_edges[order(stage, decreasing = TRUE),]
   edges <- edges[edges$values %in% label, 2:1]
   edges$values <- factor(edges$values, levels = label)

   levels(tree) <- as.character(root)
   levels(root) <- label
   levels(leaf) <- label
   levels(edges$ind) <- label

   rownames(edges) = NULL
   colnames(edges) = c("parent", "child")

   return(list(edges = edges, root = root,
               leaf = leaf, tree = tree))
}

combine_formula <- function(formula) {
   label  <- sapply(formula, attr, "label")
   nclass <- lapply(formula, attr, "nclass")
   names(nclass) <- label

   lnc <- unique(stack(nclass))
   if (!setequal(lnc$ind, label)) {
      stop("Some latent variable has not been assigned number of classes.")
   }
   lnc <- lnc[!duplicated(lnc$ind),]

   rhs <- lapply(formula, attr, "vars")
   names(rhs) <- label

   edges <- unique(utils::stack(rhs))
   edges <- edges[order(edges$ind, edges$values),]
   vars <- split(edges$values, edges$ind)
   edge <- proc_edges(edges)
   manifest <- lapply(vars[edge$leaf], function(x) x[!x %in% label])
   latent_vars <- lapply(vars, function(x) x[x %in% label])
   latent <- latent_vars[sapply(latent_vars, length) > 0]

   list(label = levels(lnc$ind),
        nclass = lnc, edge = edge,
        vars = list(manifest = manifest,
                    latent = latent))
}

constr_leaf <- function(constraints, model_table) {
   leaf <- model_table$edge$leaf
   nclass <- model_table$nclass$values
   leaf_constr <- letters[seq_len(length(leaf))]

   for (i in seq_along(constraints)) {
      constr <- constraints[[i]]
      if (!all(unlist(constr) %in% model_table$label)) next
      leaf_constr[match(constr, leaf)] <- i
   }

   nclass_leaf <- nclass[leaf[!duplicated(leaf_constr)]]

   list(leaf_constr = as.numeric(factor(leaf_constr)),
        nclass_leaf = unname(nclass_leaf),
        constr = constraints)
}

proc_formula <- function(formula, constraints) {
   formula <- sapply(formula, get_lformula)
   model_table <- combine_formula(formula)
   constr <- constr_leaf(constraints, model_table)

   list(label = model_table$label,
        nclass = model_table$nclass$values,
        root = model_table$edge$root,
        leaf = model_table$edge$leaf,
        tree = model_table$edge$tree,
        edges = model_table$edge$edges,
        vars = model_table$vars,
        constr = constr)
}
