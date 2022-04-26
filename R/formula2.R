label_formula <- function(fmla) {
   if (length(fmla) == 2) {
      return(NULL)
   }
   all.vars(fmla)[1]
}

names_latent <- function(label) {
   if (length(label) == 0) {
      return(NULL)
   }
   paste0("*", label)
}


rhs_formula <- function(fmla, labels) {
   rhs_vars  <- labels(terms(fmla))
   rhs <- paste(rhs_vars, collapse = ", ")
   rhs <- rhs_vars
   if (all(rhs_vars %in% labels)) {
      return(paste("*", rhs))
   } else if (any(rhs %in% labels)) {
      return(NULL)
   } else {
      return(rhs)
   }
}

combine_lformula <- function(formula) {
   labels <- sapply(formula, FUN = label_formula)
   rhs <- sapply(formula, FUN = rhs_formula, labels)
   names(rhs) <- labels
   assoc <- stack(rhs)


   assoc %>% group_by(ind) %>% summarise(var = paste(values, collapse = ", "))
   formula_node <- assoc[grepl("\\*", assoc$values),]
   formula_leaf <- assoc[which(!edge_index)]
}
sapply(ff, function(f) formula_label(f))

sapply(labels, function(x) rhs[[x]], simplify = FALSE)



