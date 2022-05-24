#include <Rcpp.h>
#include "auxiliary.h"
#include "transform.h"
#include "param.h"
#include "upward.h"
#include "downward.h"
#include "mstep.h"
using namespace Rcpp;

// [[Rcpp::export]]
List emFit(
   IntegerVector y, int nobs, IntegerVector nvar, List ncat,
   int nlv, int nroot, int nedge, int nleaf, int nleaf_unique,
   IntegerVector root, IntegerVector tree_index,
   IntegerVector ulv, IntegerVector vlv,
   IntegerVector leaf, IntegerVector cstr_leaf,
   IntegerVector nclass, IntegerVector nclass_leaf,
   LogicalVector init, List init_param,
   int max_iter, double tol, bool verbose, int periter = 100
) {
   int *py;

   List lst_pi(nroot), lst_tau(nedge), lst_rho(nleaf_unique);
   List lst_nrho_d(nleaf_unique), lst_nrho_n(nleaf_unique);
   std::vector<double*> ptr_pi(nroot), ptr_tau(nedge), ptr_rho(nleaf_unique);
   std::vector<double*> ptr_nrho_d(nleaf_unique), ptr_nrho_n(nleaf_unique);

   List lst_ll(nroot);
   List lst_a(nlv), lst_l(nlv), lst_j(nedge);
   std::vector<double*> ptr_ll(nroot);
   std::vector<double*> ptr_a(nlv), ptr_l(nlv), ptr_j(nedge);

   List lst_post(nlv), lst_joint(nedge);
   std::vector<double*> ptr_post(nlv), ptr_joint(nedge);

   if (init[0]) {
      List piList = init_param["pi"];
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi_init = piList[r];
         NumericVector pi = clone(pi_init);
         ptr_pi[r] = pi.begin();
         lst_pi[r] = pi;
      }
   } else {
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi = pi_gnr(nclass[root[r]]);
         ptr_pi[r] = pi.begin();
         lst_pi[r] = pi;
      }
   }

   if (init[1]) {
      List tauList = init_param["tau"];
      for (int d = 0; d < nedge; d ++) {
         NumericMatrix tau_init = tauList[d];
         NumericMatrix tau = clone(tau_init);
         ptr_tau[d] = tau.begin();
         lst_tau[d] = tau;
      }
   } else {
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d]; int v = vlv[d];
         NumericMatrix tau = tau_gnr(nclass[u], nclass[v]);
         ptr_tau[d] = tau.begin();
         lst_tau[d] = tau;
      }
   }

   if (init[2]) {
      List rhoList = init_param["rho"];
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericVector rho_init = rhoList[v];
         NumericVector lrho = clone(rho_init);
         NumericVector denom(nclass_leaf[v] * nvar[v]);
         NumericVector numer(nclass_leaf[v] * sum(ncatv));
         ptr_rho[v] = lrho.begin();
         ptr_nrho_d[v] = denom.begin();
         ptr_nrho_n[v] = numer.begin();
         lst_rho[v] = lrho;
         lst_nrho_d[v] = denom;
         lst_nrho_n[v] = numer;
      }
   } else {
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericVector lrho = rho_gnr(nclass_leaf[v], ncatv);
         NumericVector denom(nclass_leaf[v] * nvar[v]);
         NumericVector numer(nclass_leaf[v] * sum(ncatv));
         ptr_rho[v] = lrho.begin();
         ptr_nrho_d[v] = denom.begin();
         ptr_nrho_n[v] = numer.begin();
         lst_rho[v] = lrho;
         lst_nrho_d[v] = denom;
         lst_nrho_n[v] = numer;
      }
   }

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

   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;

   while ( (iter < max_iter) && (dll > tol) ) {
      iter ++;
      lastll = currll;

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


      // (maximization-step)
      // pi updates
      for (int r = 0; r < nroot; r ++)
         updatePi(ptr_pi[r], ptr_post[root[r]], nobs, nclass[root[r]]);

      // tau updates
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d]; int v = vlv[d];
         NumericMatrix ntau(nclass[u], nclass[v]);
         cumTau(ptr_joint[d], ntau.begin(), nobs, nclass[u], nclass[v]);
         updateTau(ptr_tau[d], ntau.begin(), nclass[u], nclass[v]);
      }

      // rho updates
      py = y.begin();
      for (int v = 0; v < nleaf_unique; v ++) {
         NumericVector denom = lst_nrho_d[v];
         NumericVector numer = lst_nrho_n[v];
         denom.fill(0);
         numer.fill(0);
      }
      for (int v = 0; v < nleaf; v ++) {
         int u = leaf[v];
         cumRho(ptr_nrho_d[cstr_leaf[v]], ptr_nrho_n[cstr_leaf[v]],
                py, nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]],
                nclass[u], ptr_post[u], ptr_rho[cstr_leaf[v]]);
         py += nobs * nvar[cstr_leaf[v]];
      }
      for (int v = 0; v < nleaf_unique; v ++) {
         updateRho(ptr_rho[v], ptr_nrho_n[v], ptr_nrho_d[v],
                   nobs, nclass_leaf[v], nvar[v], ncat[v]);
      }

      currll = 0;
      for (int r = 0; r < nroot; r ++) {
         double *ll = ptr_ll[r];
         for (int i = 0; i < nobs; i ++) {
            currll += ll[i];
         }
      }

      if (lastll == R_NegInf) dll = R_PosInf;
      else dll = currll - lastll;
      if (verbose) {
         Rcout << iter << " iteration  logLik " <<
            std::fixed << std::setprecision(2) << currll << "  diff ";
         if (dll > 1e-5) Rcout << std::setprecision(5) << dll;
         else Rcout << std::scientific << std::setprecision(1) << dll;
         Rcout << "      \r" << std::flush;

         if (iter % periter == 0) {
            Rcout << "\n" << std::flush;
         }
      }
   }
   if (verbose) {
      if (iter % periter != 0) Rcout << "\n";
   }


   List par;
   par["pi"] = lst_pi;
   par["tau"] = lst_tau;
   par["rho"] = lst_rho;

   List res;
   res["params"] = par;
   res["converged"] = dll < tol;
   res["niter"] = iter;

   return res;
}


// [[Rcpp::export]]
double floglik(
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

   double ll = 0;
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
      ll += calcll(ptr_l[root[r]], ptr_pi[r],
                  nobs, nclass[root[r]]);
   }

   return -ll;
}
