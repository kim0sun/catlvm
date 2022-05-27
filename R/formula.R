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

proc_edge <- function(edges) {
   label <- levels(edges$ind)
   links <- edges
   index_root <- logical(nrow(links))
   stage <- integer(nrow(links))
   while (any(!index_root)) {
      copy_links <- links
      names(copy_links) <- c("ind", "ind2")
      merged <- merge(links, copy_links, all.x = TRUE)
      rooted <- which(is.na(merged[[3]]))
      stage[rooted] <- stage[rooted] + 1
      index_root[rooted] <- TRUE
      merged$ind2[rooted] <- merged$ind[rooted]
      links <- data.frame(
         values = merged$values,
         ind = merged$ind2
      )
   }

   rt <- data.frame(root = links$ind, lv = edges$ind)
   root <- unique(rt$root)
   leaf <- unique(edges$ind[!edges$values %in% label])
   tree <- droplevels(unique(rt[order(rt$lv),])$root)
   links <- edges[order(stage, decreasing = TRUE),]
   links <- links[links$values %in% label, 2:1]
   links$values <- factor(links$values, levels = label)

   levels(tree) <- as.character(root)
   levels(root) <- label
   levels(leaf) <- label
   levels(links$ind) <- label

   rownames(links) = NULL
   colnames(links) = c("parent", "child")

   return(list(links = links, root = root,
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
   edges <- proc_edge(edges)
   manifest <- lapply(vars[edges$leaf], function(x) x[!x %in% label])
   latent_vars <- lapply(vars, function(x) x[x %in% label])
   latent <- latent_vars[sapply(latent_vars, length) > 0]

   list(label = levels(lnc$ind),
        nclass = lnc, edges = edges,
        vars = list(manifest = manifest,
                    latent = latent))
}

identify_constr <- function(constraints, model_table) {
   label <- model_table$label
   nclass <- model_table$nclass$values
   edges <- model_table$edges
   leaf <- edges$leaf
   links <- edges$links
   cstr_link <- letters[seq_len(nrow(edges$links))]
   cstr_leaf <- letters[seq_len(length(leaf))]

   for (i in seq_along(constraints)) {
      constr <- constraints[[i]]
      rm_arrow <- strsplit(constr, "->|~")
      rm_space <- lapply(rm_arrow, function(x) sub(" ", "", x))
      len <- sapply(rm_space, length)

      if (all(len == 1)) {
         ind <- sapply(rm_space, match, label)
         cstr_leaf[ind] = i
      } else if (all(len == 2)) {
         plink <- apply(edges$links, 1, paste0, collapse = "")
         ind <- sapply(rm_space, function(x)
            match(paste0(x, collapse = ""), plink))
         cstr_link[ind] = i
      } else {
         next
      }
   }

   nclass_leaf <- nclass[leaf[!duplicated(cstr_leaf)]]
   nclass_u <- nclass[links$child[!duplicated(cstr_link)]]
   nclass_v <- nclass[links$parent[!duplicated(cstr_link)]]

   list(cstr_leaf = as.numeric(factor(cstr_leaf)),
        cstr_link = as.numeric(factor(cstr_link)),
        nclass_leaf = nclass_leaf,
        nclass_u = nclass_u, nclass_v = nclass_v,
        constr = constraints)
}

proc_formula <- function(formula, constraints) {
   formula <- sapply(formula, get_lformula)
   model_table <- combine_formula(formula)
   constr <- identify_constr(constraints, model_table)

   list(label = model_table$label,
        nclass = model_table$nclass$values,
        root = model_table$edges$root,
        leaf = model_table$edges$leaf,
        tree = model_table$edges$tree,
        links = model_table$edges$links,
        vars = model_table$vars,
        constr = constr)
}
