LHSfmla <- function(lhs) {
   stopifnot(length(lhs > 1))
   list(name = lhs[[1]], nclass = lhs[[2]],
        group_var = if (length(lhs) == 3) lhs[[3]])
}

RHSfmla <- function(rhs) {
   drawItem <- function(rhs) {
      if (length(rhs) == 1) return(rhs)
      return(c(getVar(rhs[-length(rhs)][[2]]), rhs[[length(rhs)]]))
   }
   sapply(drawItem(rhs), eval, data)
}

lvFmla <- function(formula, ...) {
   if (length(formula) < 3) stop("Latent structure is wrong.")
   lc = LHSfmla(formula[[2]])
   y = model.matrix(formula[-2], data = data)[,-1]

   return(list(lc = lc[-3], y = y))
}
