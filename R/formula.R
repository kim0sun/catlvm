LHSfmla <- function(lhs) {
   stopifnot(length(lhs) > 1)
   c(name = lhs[[1]], nclass = lhs[[2]],
     group_var = if (length(lhs) == 3) lhs[[3]])
}

modelFrame <- function(...) {
   fmlrs = list(...)

   for (f in fmlrs) {
      if (length(f) < 3) stop("Latent structure is wrong.")
      lc = LHSfmla(f[[2]])
      y = all.vars(f[[3]])

      yexist = y %in% names(data)
      if (all(yexist)) msr = TRUE
      else if (all(!yexist)) msr = FALSE
      else stop("Latent structure is wrong.")
   }
}


modelFrame <- function(formula, ...) {
   if (length(formula) < 3) stop("Latent structure is wrong.")
   lc = LHSfmla(formula[[2]])
   y = all.vars(formula[[3]])

   list(lc = lc, y = y)
}


covFmla = function(formula, ...) {
   if (length(formula) < 3) stop("Model structure is wrong.")
   lc = formula[[2]][[1]]
   x = RHSfmla(formula[[3]], ...)

   list(lc = lc, x = x)
}
