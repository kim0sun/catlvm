simulate.catlvm = function(
   object, n = 500, pi = NULL, tau = NULL, rho = NULL
) {
   root = object$root
   leaf = object$leaf
   links = object$links
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
   for (d in nrow(links):1) {
      u = links[d, 1]; v = links[d, 2]
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

str_catlvm = function(formulae) {

}

