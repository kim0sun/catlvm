#include <Rcpp.h>
#include "upward.h"
#include "downward.h"
#include "transform.h"
using namespace Rcpp;

// [[Rcpp::export]]
List calcModel(
      List param, IntegerVector y, int nobs,
      IntegerVector nvar, List ncat,
      int nlv, int nroot, int nlink, int nleaf,
      int nlink_unique, int nleaf_unique,
      IntegerVector root, IntegerVector tree_index,
      IntegerVector ulv, IntegerVector vlv, IntegerVector leaf,
      IntegerVector cstr_link, IntegerVector cstr_leaf,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v
) {
   List lst_pi = param["pi"];
   List lst_tau = param["tau"];
   List lst_rho = param["rho"];
   std::vector<double*> ptr_pi(nroot);
   std::vector<double*> ptr_tau(nlink_unique);
   std::vector<double*> ptr_rho(nleaf_unique);

   List score_pi(nroot);
   List score_tau(nlink_unique);
   List score_rho(nleaf_unique);
   std::vector<double*> ptr_spi(nroot);
   std::vector<double*> ptr_stau(nlink_unique);
   std::vector<double*> ptr_srho(nleaf_unique);

   for (int r = 0; r < nroot; r ++) {
      int nk = nclass[root[r]];
      NumericVector pi = lst_pi[r];
      NumericMatrix s_pi(nk - 1, nobs);
      ptr_pi[r] = pi.begin();
      score_pi[r] = s_pi;
      ptr_spi[r] = s_pi.begin();
   }

   for (int d = 0; d < nlink_unique; d ++) {
      int nk = nclass_u[d], nl = nclass_v[d];
      NumericMatrix tau = lst_tau[d];
      NumericMatrix s_tau(nl * (nk - 1), nobs);
      ptr_tau[d] = tau.begin();
      score_tau[d] = s_tau;
      ptr_stau[d] = s_tau.begin();
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      int nk = nclass_leaf[v];
      IntegerVector ncatv = ncat[v];
      NumericVector rho = lst_rho[v];
      NumericMatrix s_rho(nk * sum(ncatv - 1), nobs);
      ptr_rho[v] = rho.begin();
      score_rho[v] = s_rho;
      ptr_srho[v] = s_rho.begin();
   }

   NumericVector ll(nobs);
   List lst_ll(nroot);
   List lst_a(nlv), lst_l(nlv), lst_j(nlink);
   std::vector<double*> ptr_ll(nroot);
   std::vector<double*> ptr_a(nlv), ptr_l(nlv), ptr_j(nlink);

   List lst_post(nlv), lst_joint(nlink);
   std::vector<double*> ptr_post(nlv), ptr_joint(nlink);

   for (int r = 0; r < nroot; r ++) {
      NumericVector ll(nobs);
      lst_ll[r] = ll;
      ptr_ll[r] = ll.begin();
   }

   for (int d = 0; d < nlink; d ++) {
      int u = ulv[d], v = vlv[d];
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

   // (expectation-step)
   // initiate lambda
   int *y_ = y.begin();
   for (int v = 0; v < nleaf; v ++) {
      upInit(y_, ptr_rho[cstr_leaf[v]],
             ptr_l[leaf[v]], nclass[leaf[v]],
             nobs, nvar[cstr_leaf[v]],
             ncat[cstr_leaf[v]]);
      y_ += nobs * nvar[cstr_leaf[v]];
   }

   // upward recursion
   for (int d = nlink - 1; d > -1; d --) {
      int u = ulv[d], v = vlv[d];
      upRec(ptr_l[v], ptr_j[d], ptr_l[u],
            ptr_tau[cstr_link[d]],
            nobs, nclass[u], nclass[v]);
   }

   // calculate loglik
   for (int r = 0; r < nroot; r ++) {
      calclli(ptr_l[root[r]], ptr_pi[r], ll.begin(),
              nobs, nclass[root[r]]);
   }

   // initiate alpha
   for (int r = 0; r < nroot; r ++) {
      dnInit(ptr_a[root[r]], ptr_l[root[r]],
             ptr_pi[r], ptr_post[root[r]],
             ptr_ll[r], nobs, nclass[root[r]]);
   }

   // Downward recursion
   for (int d = 0; d < nlink; d ++) {
      int u = ulv[d], v = vlv[d];
      dnRec(ptr_a[u], ptr_a[v],
            ptr_l[u], ptr_l[v], ptr_j[d],
            nobs, nclass[u], nclass[v],
            ptr_tau[cstr_link[d]],
            ptr_post[u], ptr_joint[d],
            ptr_ll[tree_index[d]]);
   }


   for (int r = 0; r < nroot; r ++) {
      int nk = nclass[root[r]];
      double *score_ = ptr_spi[r];
      double *post_ = ptr_post[root[r]];
      double *pi_ = ptr_pi[r];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < nk - 1; k ++)
            score_[k] = exp(post_[k]) - exp(pi_[k]);
         score_ += nk - 1;
         post_  += nk;
      }
   }

   for (int d = 0; d < nlink; d ++) {
      int u = ulv[d], v = vlv[d];
      double *score_ = ptr_stau[cstr_link[d]];
      double *joint_ = ptr_joint[d];
      double *post_ = ptr_post[v];
      for (int i = 0; i < nobs; i ++) {
         double *tau_ = ptr_tau[cstr_link[d]];
         for (int l = 0; l < nclass[v]; l ++) {
            for (int k = 0; k < nclass[u] - 1; k ++)
               score_[k] += exp(joint_[k]) - exp(tau_[k] + post_[l]);
            score_ += nclass[u] - 1;
            joint_ += nclass[u];
            tau_ += nclass[u];
         }
         post_ += nclass[v];
      }
   }

   for (int v = 0; v < nleaf; v ++) {
      int nk = nclass[leaf[v]];
      IntegerVector ncatv = ncat[cstr_leaf[v]];
      double *score_ = ptr_srho[cstr_leaf[v]];
      double *post_ = ptr_post[v];
      for (int i = 0; i < nobs; i ++) {
         double *rho_ = ptr_rho[cstr_leaf[v]];
         for (int k = 0; k < nk; k ++) {
            for (int m = 0; m < nvar[cstr_leaf[v]]; m ++) {
               // if (y_[m] < ncatv[m] - 1)
               //    score_[y_[m] - 1] += exp(post_[k]);
               for (int r = 0; r < ncatv[m] - 1; r ++)
                  score_[r] -= exp(rho_[r] + post_[k]);
               score_ += ncatv[m] - 1;
               rho_ += ncatv[m];
            }
         }
         post_ += nk;
         y_ += nvar[cstr_leaf[v]];
      }
   }

   // Posterior transpose
   for (int v = 0; v < nlv; v ++) {
      NumericMatrix post = lst_post[v];
      lst_post[v] = transpose(post);
   }
   for (int d = 0; d < nlink; d ++) {
      NumericMatrix joint = lst_joint[d];
      lst_joint[d] = transpose(joint);
   }

   List score;
   score["pi"] = score_pi;
   score["tau"] = score_tau;
   score["rho"] = score_rho;

   List res;
   res["ll"] = ll;
   res["score"] = score;
   res["post"] = lst_post;
   res["joint"] = lst_joint;
   res["lambda"] = lst_l;

   return res;
}
