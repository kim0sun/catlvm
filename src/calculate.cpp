#include <Rcpp.h>
#include "upward.h"
#include "downward.h"
#include "transform.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calcll(
      NumericVector param,
      IntegerVector y, int nobs, IntegerVector nvar, List ncat,
      int nlv, int nroot, int nedge, int nleaf, int nleaf_unique,
      IntegerVector root, IntegerVector ulv, IntegerVector vlv,
      IntegerVector leaf, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf
) {
   int *py;
   List lst_pi(nroot);
   List lst_tau(nedge);
   List lst_rho(nleaf_unique);
   std::vector<double*> ptr_pi(nroot);
   std::vector<double*> ptr_tau(nedge);
   std::vector<double*> ptr_rho(nleaf_unique);

   List lst_l(nlv);
   std::vector<double*> ptr_l(nlv);
   double *param_ = param.begin();

   NumericVector ll(nobs);

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi = logistic_pi(param_, nclass[root[r]]);
      ptr_pi[r] = pi.begin();
      lst_pi[r] = pi;
      param_ += nclass[root[r]] - 1;
   }

   for (int d = 0; d < nedge; d ++) {
      int u = ulv[d]; int v = vlv[d];
      NumericMatrix tau = logistic_tau(param_, nclass[u], nclass[v]);
      ptr_tau[d] = tau.begin();
      lst_tau[d] = tau;
      param_ += nclass[v] * (nclass[u] - 1);
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho = logistic_rho(param_, nclass_leaf[v], ncatv);
      ptr_rho[v] = rho.begin();
      lst_rho[v] = rho;
      param_ += nclass_leaf[v] * sum(ncatv - 1);
   }

   for (int v = 0; v < nlv; v ++) {
      NumericMatrix lambda(nclass[v], nobs);
      ptr_l[v] = lambda.begin();
      lst_l[v] = lambda;
   }

   // initiate lambda
   py = y.begin();
   for (int v = 0; v < nleaf; v ++) {
      upInit(py, ptr_rho[cstr_leaf[v]], ptr_l[leaf[v]],
             nclass[leaf[v]], nobs,
             nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
      py += nobs * nvar[cstr_leaf[v]];
   }

   // upward recursion
   for (int d = nedge - 1; d > -1; d --) {
      int u = ulv[d];
      int v = vlv[d];
      upRec2(ptr_l[v], ptr_l[u], ptr_tau[d],
             nobs, nclass[u], nclass[v]);
   }

   for (int r = 0; r < nroot; r ++) {
      calclli(ptr_l[root[r]], ptr_pi[r], ll.begin(),
              nobs, nclass[root[r]]);
   }

   return ll;
}


// [[Rcpp::export]]
List calcModel(
      List param, IntegerVector y, int nobs, IntegerVector nvar, List ncat,
      int nlv, int nroot, int nedge, int nleaf, int nleaf_unique,
      IntegerVector root, IntegerVector tree_index,
      IntegerVector ulv, IntegerVector vlv,
      IntegerVector leaf, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf
) {
   int *py;

   List lst_pi = param["pi"], lst_tau = param["tau"], lst_rho = param["rho"];
   std::vector<double*> ptr_pi(nroot), ptr_tau(nedge), ptr_rho(nleaf_unique);

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi = lst_pi[r];
      ptr_pi[r] = pi.begin();
   }

   for (int d = 0; d < nedge; d ++) {
      NumericMatrix tau = lst_tau[d];
      ptr_tau[d] = tau.begin();
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      NumericVector rho = lst_rho[v];
      ptr_rho[v] = rho.begin();
   }

   NumericVector lls(nobs);
   List lst_ll(nroot);
   List lst_a(nlv), lst_l(nlv), lst_j(nedge);
   std::vector<double*> ptr_ll(nroot);
   std::vector<double*> ptr_a(nlv), ptr_l(nlv), ptr_j(nedge);

   List lst_post(nlv), lst_joint(nedge);
   std::vector<double*> ptr_post(nlv), ptr_joint(nedge);

   for (int r = 0; r < nroot; r ++) {
      NumericVector ll(nobs);
      lst_ll[r] = ll;
      ptr_ll[r] = ll.begin();
   }

   for (int d = 0; d < nedge; d ++) {
      int u = ulv[d]; int v = vlv[d];
      NumericMatrix joint(nclass[u] * nclass[v], nobs);
      ptr_joint[d] = joint.begin();
      lst_joint[d] = joint;
      NumericMatrix jlambda(nclass[v], nobs);
      ptr_j[d] = jlambda.begin();
      lst_j[d] = jlambda;
   }

   for (int v = 0; v < nlv; v ++) {
      NumericMatrix post(nclass[v], nobs);
      NumericMatrix alpha(nclass[v], nobs);
      NumericMatrix lambda(nclass[v], nobs);
      ptr_post[v] = post.begin();
      ptr_a[v] = alpha.begin();
      ptr_l[v] = lambda.begin();
      lst_post[v] = post;
      lst_a[v] = alpha;
      lst_l[v] = lambda;
   }

   // lambda, cleaning
   for (int v = 0; v < nlv; v ++) {
      NumericMatrix lambda = lst_l[v];
      lambda.fill(0);
   }

   // (expectation-step)
   // initiate lambda
   py = y.begin();
   for (int v = 0; v < nleaf; v ++) {
      upInit(py, ptr_rho[cstr_leaf[v]], ptr_l[leaf[v]], nclass[leaf[v]],
             nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
      py += nobs * nvar[cstr_leaf[v]];
   }

   // upward recursion
   for (int d = nedge - 1; d > -1; d --) {
      int u = ulv[d];
      int v = vlv[d];
      upRec(ptr_l[v], ptr_j[d], ptr_l[u], ptr_tau[d],
            nobs, nclass[u], nclass[v]);
   }

   // calculate loglik
   for (int r = 0; r < nroot; r ++) {
      calclli(ptr_l[root[r]], ptr_pi[r], lls.begin(),
              nobs, nclass[root[r]]);
   }

   // initiate alpha
   for (int r = 0; r < nroot; r ++) {
      dnInit(ptr_a[root[r]], ptr_l[root[r]],
             ptr_pi[r], ptr_post[root[r]],
             ptr_ll[r], nobs, nclass[root[r]]);
   }

   // Downward recursion
   for (int d = 0; d < nedge; d ++) {
      int u = ulv[d];
      int v = vlv[d];
      dnRec(ptr_a[u], ptr_a[v], ptr_l[u], ptr_l[v], ptr_j[d],
            nobs, nclass[u], nclass[v], ptr_tau[d], ptr_post[u],
            ptr_joint[d], ptr_ll[tree_index[d]]);
   }


   List res;
   res["ll"] = lls;
   res["post"] = lst_post;
   res["lambda"] = lst_l;

   return res;
}
