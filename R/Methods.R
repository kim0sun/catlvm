simulate.catlvm = function(
   object, n = 500, pi = NULL, tau = NULL, rho = NULL
) {
   root = object$root
   leaf = object$leaf
   ulv = object$ulv
   vlv = object$vlv
   cstr_lf = object$cstr_lf
   cstr_lk = object$cstr_lk

   nc = object$nc
   unc = object$unc
   lnc = object$lnc
   nclf = object$nclf

   nleaf = object$nleaf
   nuleaf = object$nuleaf
   nlink = object$nlink
   nulink = object$nulink

   nvar = object$nvar
   ncat = object$ncat

   if (is.null(pi)) pi = exp(pi_gnr(nc[root], 1))
   if (is.null(tau)) {
      tau = list()
      for (d in seq(nulink)) {
         tau[[d]] = exp(matrix(tau_gnr(unc[d], lnc[d], 1), unc[d]))
      }
   }
   if (is.null(rho)) {
      rho = list()
      for (v in seq(nuleaf)) {
         rho[[v]] = exp(rho_gnr(nclf[v], ncat[[v]]))
      }
   }

   cls = list()
   cls[[root]] = root_gnr(n, nc[root], pi)
   for (d in rev(seq(ulv))) {
      u = ulv[d]; v = vlv[d]
      tau_d = tau[[cstr_lk[d]]]
      cls[[u]] = cls_gnr(n, nc[u], nc[v], cls[[v]], tau_d);
   }
   y = list()
   for (v in seq(nleaf)) {
      rho_v = rho[[cstr_lf[v]]]
      y[[v]] = y_gnr(n, nclf[cstr_lf[v]], ncat[[cstr_lf[v]]], cls[[leaf[v]]], rho_v)
   }

   y
}

catlvm = function(
   measurementModel = NULL,
   latentStructure = NULL,
   constraints = NULL,
   data = NULL
) {

   class(model) = "catlvm"
}

estimate.catlvm = function(
   x,  ...
) {

   class(res) = c("catlvm", "catlvmFitted")
}

npar.catlvm = function(object, ...) {

}

nobs.catlvm = function(object, ...) {

}

print.catlvm = function(x, ...) {

}

summary.catlvm = function(object, ...) {

}

logLik.catlvm = function(object, ...) {

}

coef.catlvm = function(object, ...) {

}

vcov.catlvm = function(object, ...) {

}
