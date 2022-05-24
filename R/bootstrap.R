bootstrap.dcatlvm <- function(
   param, method, args, B = 100,
   max.iter = 500, verbose = FALSE
) {
   for (b in seq_len(B)) {
      ysim <- ysim(args$nobs, args$nlv, args$root - 1, args$leaf - 1,
                   args$u - 1, args$v - 1, args$nclass, args$nroot,
                   args$nleaf, args$nedge, args$ncat, args$cstr_leaf - 1,
                   param$pi, param$tau, param$rho, FALSE)$y
      lparam <- logit_param(param, args)
      nlm_fit <- nlm(floglik, unlist(lparam),
         y = ysim, nobs = args$nobs, nvar = args$nvar, ncat = args$ncat,
         nlv = args$nlv, nroot = args$nroot, nedge = args$nedge,
         nleaf = args$nleaf, nleaf_unique = args$nleaf_unique,
         root = args$root - 1, ulv = args$u - 1, vlv = args$v - 1,
         leaf = args$leaf - 1, cstr_leaf = args$cstr_leaf - 1,
         nclass = args$nclass, nclass_leaf = args$nclass_leaf,
         iterlim = max.iter, steptol = 1e-6
      )

   }
}
