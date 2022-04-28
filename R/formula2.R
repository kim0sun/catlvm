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
dput(fmla[[2]][[1]])
f = L1|3 ~ X1 + X2
f[[2]][[3]]
setdiff(as.numeric(fmla[[2]]), label_formula(fmla))

label_formula(f)
f = L1 + L2 ~ X1 + X2
length(f[[2]])

fmla = L1[3] ~ X1 + X2
fmla[[2]][[2]]

attributes(fmla)

rhs_formula <- function(fmla, labels) {

}

get_lformula <- function(fmla) {
   if (!inherits(fmla, "formula"))
      fmla <- as.formula(fmla)

   attr(fmla, "label") = label_formula(fmla)
   attr(fmla, "nclass") = nclass_formula(fmla)
   attr(fmla, "vars") = labels(terms(fmla))

   fmla
}
formula = list(
   L1[3] ~ X1 + X2,
   L1 + L2 ~ X1,
   L2[3] ~ X3,
   L2[3] ~ X4,
   P ~ L1 + L2,
   PF[3] ~ L1 + L2
)

combine_formula <- function(formula) {
   formula <- sapply(formula, get_lformula)
   rhs <- lapply(formula, attr, "vars")
   names(rhs) <- sapply(formula, attr, "label")
}

proc_formula <- function(formula) {
   formula <- sapply(formula, get_lformula)
   label <- sapply(formula, label_formula)

   formula_reg <- list()
   formula_str <- list()
   for (f in formula) {
      rl <- all(attr(fmla, "vars") %in% labels)
      n0 <- length(attr(fmla, "nclass")) == 0
      if (ro && no)
         formula_reg <- append(formula_reg, f)
      if ()
         formula_str <- append(formula_leaf, f)
   }
   rhs <- lapply(formula_leaf)
   names(rhs) <- labels
   assoc <- stack(rhs)
}

regression_lhs <- function(fmla) {
   if (length(fmla) == 2) {
      return(FALSE)
   } else if (length(fmla[[2]]) == 1) {
      return(TRUE)
   } else {
      reutrn(FALSE)
   }
}

classify_formula <- function(formula) {
   labels <- sapply(formula, FUN = label_formula)
   regress <-
   latent  <-
}

combine_lformula <- function(formula) {
   labels <- sapply(formula, FUN = label_formula)
   rhs <- sapply(formula, FUN = rhs_formula, labels)
   names(rhs) <- labels
   assoc <- stack(rhs)

   formula_node <- assoc[grepl("\\*", assoc$values),]
   formula_leaf <- assoc[which(!edge_index)]
}
sapply(ff, function(f) formula_label(f))

sapply(labels, function(x) rhs[[x]], simplify = FALSE)



