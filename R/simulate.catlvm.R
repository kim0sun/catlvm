simulate.catlvm = function(
   object, n = 500, pi, tau, rho
) {
   root = struct$root
   leaf = struct$leaf
   links = struct$links
   cstr_lf = struct$cstr_lf
   cstr_lk = struct$cstr_lk

   nc = struct$nc
   unc = struct$unc
   lnc = struct$lnc
   nclf = struct$nclf

   nleaf = struct$nleaf
   nuleaf = struct$nuleaf
   nlink = struct$nlink
   nulink = struct$nulink

   nvar = struct$nvar
   ncat = struct_ncat

   if (is.null(pi)) pi = pi_gnr(nc[root], 1)
   if (is.null(tau)) {
      tau = list()
      for (d in seq(nulink)) {
         tau[[d]] = matrix(tau_gnr(unc[d], lnc[d], 1), unc[d])
      }
   }
   if (is.null(rho)) {
      rho = list()
      for (v in seq(nuleaf)) {
         rho[[v]] = rho_gnr(nclf[v], ncat[v])
      }
   }

   cls = list()
   cls[[root]] = root_gnr(n, nc[root], pi)
   for (d in nrow(links):1) {
      u = links[d, 1]; v = links[d, 2]
      tau_d = tau[[cstr_lk[d]]]
      cls[[u]] = cls_gnr(n, nc[u], nc[v], tau_d);
   }
   y = list()
   for (v in seq(nleaf)) {
      rho_v = rho[[cstr_lf[v]]]
      y[[leaf[v]]] = y_gnr(n, nc[leaf[v]], ncat, cls[[leaf[v]]], rho_v)
   }
}

str_catlvm = function(formulae) {

}

