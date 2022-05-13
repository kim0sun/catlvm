#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

int sample1(int n, double *prob) {
   double ran = R::runif(0, 1);
   double cp = 0;
   for (int i = 0; i < n; i ++) {
      cp += exp(prob[i]);
      if (ran < cp) return i;
   }
   return n - 1;
}

// [[Rcpp::export]]
NumericVector logit_pi(
   NumericVector pi, int nclass
) {
   NumericVector lpi(nclass - 1);

   for (int k = 0; k < nclass - 1; k ++)
      lpi[k] = pi[k] - pi[nclass - 1];

   return lpi;
}

// [[Rcpp::export]]
NumericMatrix logit_tau(
   NumericMatrix tau, int nk, int nl
) {
   NumericMatrix ltau(nk - 1, nl);
   double *tau_ = tau.begin();
   double *ltau_ = ltau.begin();

   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk - 1; k ++) {
         ltau_[k] = tau_[k] - tau_[nk - 1];
      }
      ltau_ += nk - 1;
      tau_  += nk;
   }

   return ltau;
}


// [[Rcpp::export]]
NumericVector logit_rho(
   NumericVector rho, IntegerVector ncat, int nclass
) {
   NumericVector lrho(nclass * sum(ncat - 1));
   double *rho_ = rho.begin();
   double *lrho_ = lrho.begin();

   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         for (int r = 0; r < ncat[m] - 1; r ++) {
            lrho_[r] = rho_[r] - rho_[ncat[m] - 1];
         }
         lrho_ += ncat[m] - 1;
         rho_  += ncat[m];
      }
   }
   return lrho;
}


NumericVector logistic_pi(double *lpi, int nclass) {
   NumericVector pi(nclass);
   double denom = 1;
   for (int i = 0; i < nclass - 1; i ++) {
      pi[i] = lpi[i];
      denom += exp(lpi[i]);
   }
   for (int i = 0; i < nclass; i ++)
      pi[i] -= log(denom);

   return pi;
}

NumericMatrix logistic_tau(double *ltau, int nk, int nl) {
   NumericMatrix tau(nk, nl);
   double *tau_ = tau.begin();
   double denom = 0;
   for (int l = 0; l < nl; l ++) {
      denom = 1;
      for (int k = 0; k < nk - 1; k ++) {
         tau_[k] = ltau[k];
         denom += exp(ltau[k]);
      }
      for (int k = 0; k < nk; k ++) {
         tau_[k] -= log(denom);
      }
      tau_ += nk;
      ltau += nk - 1;
   }

   return tau;
}

NumericVector logistic_rho(double *lrho, int nclass, IntegerVector ncat) {
   NumericVector rho(nclass * sum(ncat));
   double *rho_ = rho.begin();
   double denom = 0;
   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         denom = 1;
         for (int r = 0; r < ncat[m] - 1; r ++) {
            rho_[r] = lrho[r];
            denom += exp(lrho[r]);
         }
         for (int r = 0; r < ncat[m]; r ++) {
            rho_[r] -= log(denom);
         }
         rho_ += ncat[m];
         lrho += ncat[m] - 1;
      }
   }

   return rho;
}


// [[Rcpp::export]]
NumericVector elogdiri(NumericVector a) {
   return digamma(a) - R::digamma(sum(a));
}

// [[Rcpp::export]]
NumericVector plogdiri(
   NumericVector a, NumericVector b
) {
   return lgammaf(sum(a)) - lgamma(a) +
      sum( (a - 1) * elogdiri(b) );;
}

// [[Rcpp::export]]
NumericVector pi_gnr(int nk) {
   NumericVector pi(nk);
   pi.fill(-log(nk));
   return pi;
}

// [[Rcpp::export]]
NumericMatrix tau_gnr(int nk, int nl) {
   NumericMatrix tau(nk, nl);
   double *ptau = tau.begin();

   for (int l = 0; l < nl; l ++) {
      NumericVector p = runif(nk, 0, 1);
      for (int k = 0; k < nk; k ++) {
         ptau[k] = log(p[k]) - log(sum(p));
      }
      ptau += nk;
   }

   return tau;
}

void prev_gnr(
   double *prev, double *upr, double *tau,
   int nobs, int nk, int nl
) {
   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk; k ++)
         prev[k] += exp(tau[k] + upr[l]);
      tau += nk;
   }
}


// [[Rcpp::export]]
NumericVector rho_gnr(int nk, IntegerVector ncat) {
   NumericVector pv(nk * sum(ncat));
   double *ppv = pv.begin();

   for (int k = 0; k < nk; k ++) {
      for (int m = 0; m < ncat.length(); m ++) {
         NumericVector p = runif(ncat[m], 0, 1);
         for (int r = 0; r < ncat[m]; r ++) {
            ppv[r] = log(p[r] / sum(p));
         }
         ppv += ncat[m];
      }
   }
   return pv;
}

// [[Rcpp::export]]
IntegerVector root_gnr(
   int nobs, int nk,
   Nullable<NumericVector> prob = R_NilValue
) {
   NumericVector pi;
   if (prob.isNull()) {
      pi = pi_gnr(nk);
   } else {
      pi = as<NumericVector>(prob);
   }

   IntegerVector cls(nobs);
   for (int i = 0; i < nobs; i ++) {
      cls[i] = sample1(nk, pi.begin());
   }

   return cls;
}

// [[Rcpp::export]]
IntegerVector cls_gnr(
   int nobs, int nk, int nl, IntegerVector v,
   Nullable<NumericMatrix> prob = R_NilValue
) {
   NumericMatrix tau;
   if (prob.isNull()) {
      tau = tau_gnr(nk, nl);
   } else {
      tau = as<NumericMatrix>(prob);
   }

   IntegerVector cls(nobs);
   double *ptau = tau.begin();
   for (int i = 0; i < nobs; i ++) {
      cls[i] = sample1(nk, ptau + v[i] * nk);
   }

   return cls;
}

// [[Rcpp::export]]
IntegerMatrix y_gnr(
   int nobs, int nk, IntegerVector ncat,
   IntegerVector cls,
   Nullable<NumericVector> prob = R_NilValue
) {
   NumericVector rho;
   if (prob.isNull()) {
      rho = rho_gnr(nk, ncat);
   } else {
      rho = as<NumericVector>(prob);
   }

   int nvar = ncat.length();
   IntegerMatrix y(nvar, nobs);
   for (int i = 0; i < nobs; i ++) {
      double *pos = rho.begin() + cls[i] * sum(ncat);
      for (int m = 0; m < nvar; m ++) {
         y[i * nvar + m] = sample1(ncat[m], pos) + 1;
         pos += ncat[m];
      }
   }

   return y;
}

// [[Rcpp::export]]
List ysim(
   int nsim, int nlv, IntegerVector root, IntegerVector leaf,
   Nullable<IntegerVector> ulv, Nullable<IntegerVector> vlv,
   IntegerVector nclass, int nroot, int nleaf, int nedge, List ncat,
   IntegerVector cstr_leaf, List pi, List tau, List rho, bool print_class
) {
   List cls(nlv);

   for (int r = 0; r < nroot; r ++)
      cls[root[r]] = root_gnr(nsim, nclass[root[r]], pi[r]);

   if (ulv.isNotNull() && vlv.isNotNull()) {
      IntegerVector u = as<IntegerVector>(ulv);
      IntegerVector v = as<IntegerVector>(vlv);
      for (int d = 0; d < nedge; d ++) {
         cls[u[d]] = cls_gnr(nsim, nclass[u[d]], nclass[v[d]],
                             cls[v[d]], tau[d]);
      }
   }

   List y(nleaf);
   for (int v = 0; v < nleaf; v ++)
      y[v] = y_gnr(nsim, nclass[leaf[v]], ncat[cstr_leaf[v]],
                   cls[leaf[v]], rho[cstr_leaf[v]]);

   if (print_class) {
      List ret;
      ret["y"] = y;
      ret["class"] = cls;
      return ret;
   }

   return y;
}

// void calcSat(
//    int *ymiss, int nobs, int nvar, IntegerVector ncat,
//    double *ccount
// ) {
//    int npattern = 1;
//    IntegerVector jumps(nvar);
//    for (int m = 0; m < nvar; m ++) {
//       jumps[m] = npattern;
//       npattern *= ncat[m];
//    }
//
//    int ind = 0;
//    for (int m = 0; m < nvar; m ++) {
//       if (y[m] > 0) {
//          for (int r = 0; r < ncat[m]; r ++)
//             ind += jumps[m] * (y[m] - 1);
//       } else {
//
//       }
//    }
//    ccount[ind] += 1;
//
//    ccount.fill(0);
// }

// [[Rcpp::export]]
double calcSll(
   double *cellp, int *cellc, int ncell
) {
   double ll;
   for (int i = 0; i < ncell; i ++) {
      ll += cellc[i] * log(cellp[i]);
   }
   return ll;
}

// Upward-Recursion
void upInit(
   int *y, const double *ptr_rho, double *lambda,
   int nk, int nobs, int nvar, IntegerVector ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *rho = (double *) ptr_rho;
      for (int k = 0; k < nk; k ++) {
         for (int m = 0; m < nvar; m ++) {
            if (y[m] > 0)
               lambda[k] += rho[y[m] - 1];
            rho += ncat[m];
         }
      }
      y += nvar;
      lambda += nk;
   }
}

void upRec(
   double *lambda, double *jlambda, double *llambda,
   const double *tau, int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *tau_ = (double*) tau;
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(tau_[k] + llambda[k]);

         tau_ += nk;
         jlambda[l] = log(ml);
         lambda[l] += log(ml);
      }
      llambda += nk;
      jlambda += nl;
      lambda  += nl;
   }
}

void upRec2(
   double *lambda, double *llambda, const double *tau,
   int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *tau_ = (double*) tau;
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(tau_[k] + llambda[k]);

         tau_ += nk;
         lambda[l] += log(ml);
      }
      llambda += nk;
      lambda  += nl;
   }
}

// Downward-Recursion
void dnInit(
   double *alpha, double *lambda, double *pi,
   double *post, double *ll, int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         alpha[k] = pi[k];
         post[k] = alpha[k] + lambda[k];
         lik += exp(post[k]);
      }
      ll[i] = log(lik);

      for (int k = 0; k < nclass; k ++) {
         post[k] -= ll[i];
      }

      alpha += nclass;
      lambda += nclass;
      post += nclass;
   }
}

double calcll(
   double *lambda, double *pi,
   int nobs, int nclass
) {
   double ll = 0;
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         lik += exp(pi[k] + lambda[k]);
      }
      ll += log(lik);

      lambda += nclass;
   }
   return ll;
}

void calclli(
   double *lambda, double *pi,
   double *ll, int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         lik += exp(pi[k] + lambda[k]);
      }
      ll[i] += log(lik);

      lambda += nclass;
   }
}

void dnRec(
   double *alpha, double *ualpha,
   double *lambda, double *ulambda, double *jlambda,
   int nobs, int nk, int nl, double *tau,
   double *post, double *joint, double *ll
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         double svl = 0;
         for (int l = 0; l < nl; l ++) {
            val = tau[k + l * nk] + ualpha[l] + ulambda[l] - jlambda[l];
            joint[k + l * nk] = val + lambda[k] - ll[i];
            svl += exp(val);
         }
         alpha[k] = log(svl);
         post[k] = alpha[k] + lambda[k] - ll[i];
      }

      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      lambda += nk; ulambda += nl; jlambda += nl;
   }
}

// M-step
void updatePi(
   double *pi, double *post, int nobs, int nclass
) {
   NumericVector npi(nclass);
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++) {
         npi[k] += exp(post[k]);
      }
      post += nclass;
   }
   for (int k = 0; k < nclass; k ++) {
      pi[k] = log(npi[k]) - log(sum(npi));
   }
}

void cumTau(
   double *joint, const double *ntau,
   int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *ntau_ = (double *) ntau;
      for (int l = 0; l < nl; l ++) {
         for (int k = 0; k < nk; k ++) {
            ntau_[k] += exp(joint[k]);
         }
         joint += nk;
         ntau_ += nk;
      }
   }
}

void updateTau(
   double *tau, double *ntau,
   int nk, int nl
) {
   for (int l = 0; l < nl; l ++) {
      double stau = 0;
      for (int k = 0; k < nk; k ++) {
         stau += ntau[k];
      }
      for (int k = 0; k < nk; k ++) {
         tau[k] = log(ntau[k]) - log(stau);
      }
      ntau += nk;
      tau  += nk;
   }
}


void cumRho(
   const double *denom, const double *numer,
   int *y, int nobs, int nvar, IntegerVector ncat,
   int nk, double *post, const double *old_rho
) {
   for (int i = 0; i < nobs; i ++) {
      double *dnm = (double *) denom;
      double *nmr = (double *) numer;
      double *rho = (double *) old_rho;
      for (int k = 0; k < nk; k ++) {
         for (int m = 0; m < nvar; m ++) {
            dnm[m] += exp(post[k]);
            if (y[m] > 0) nmr[y[m] - 1] += exp(post[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  nmr[r] += exp(post[k] + rho[r]);
               }
            }
            nmr += ncat[m];
            rho += ncat[m];
         }
         dnm += nvar;
      }
      post += nk;
      y += nvar;
   }
}

void updateRho(
   double *rho, double *numer, double *denom,
   int nobs, int nclass, int nvar, IntegerVector ncat
) {
   for (int k = 0; k < nclass; k ++) {
      for (int m = 0; m < nvar; m ++) {
         for (int r = 0; r < ncat[m]; r ++) {
            rho[r] = log(numer[r]) - log(denom[m]);
         }
         rho   += ncat[m];
         numer += ncat[m];
      }
      denom += nvar;
   }
}


// [[Rcpp::export]]
List par_gnr(
   int nobs, IntegerVector nvar, List ncat,
   int nroot, int nedge, int nleaf_unique,
   IntegerVector root, IntegerVector ulv, IntegerVector vlv,
   IntegerVector nclass, IntegerVector nclass_leaf,
   LogicalVector init, List init_param
) {
   List lst_pi(nroot);
   List lst_tau(nedge);
   List lst_rho(nleaf_unique);

   if (init[0]) {
      List piList = init_param["pi"];
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi_init = piList[r];
         NumericVector pi = clone(pi_init);
         lst_pi[r] = pi;
      }
   } else {
      for (int r = 0; r < nroot; r ++) {
         NumericVector pi = pi_gnr(nclass[root[r]]);
         lst_pi[r] = pi;
      }
   }

   if (init[1]) {
      List tauList = init_param["tau"];
      for (int d = 0; d < nedge; d ++) {
         NumericMatrix tau_init = tauList[d];
         NumericMatrix tau = clone(tau_init);
         lst_tau[d] = tau;
      }
   } else {
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d]; int v = vlv[d];
         NumericMatrix tau = tau_gnr(nclass[u], nclass[v]);
         lst_tau[d] = tau;
      }
   }

   if (init[2]) {
      List rhoList = init_param["rho"];
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericVector rho_init = rhoList[v];
         NumericVector lrho = clone(rho_init);
         lst_rho[v] = lrho;
      }
   } else {
      for (int v = 0; v < nleaf_unique; v ++) {
         IntegerVector ncatv = ncat[v];
         NumericVector lrho = rho_gnr(nclass_leaf[v], ncatv);
         lst_rho[v] = lrho;
      }
   }

   List par;
   par["pi"] = lst_pi;
   par["tau"] = lst_tau;
   par["rho"] = lst_rho;

   return par;
}

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


// [[Rcpp::export]]
NumericVector calclli(
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
List calcPost(
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

   return lst_post;
}


// [[Rcpp::export]]
List logit2log(
   NumericVector param, int nobs, List ncat,
   int nroot, int nedge, int nleaf_unique,
   IntegerVector root, IntegerVector ulv, IntegerVector vlv,
   IntegerVector nclass, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nedge);
   List lst_rho(nleaf_unique);

   double *param_ = param.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi = logistic_pi(param_, nclass[root[r]]);
      lst_pi[r] = pi;
      param_ += nclass[root[r]] - 1;
   }

   for (int d = 0; d < nedge; d ++) {
      int u = ulv[d]; int v = vlv[d];
      NumericMatrix tau = logistic_tau(param_, nclass[u], nclass[v]);
      lst_tau[d] = tau;
      param_ += nclass[v] * (nclass[u] - 1);
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho = logistic_rho(param_, nclass_leaf[v], ncatv);
      lst_rho[v] = rho;
      param_ += nclass_leaf[v] * sum(ncatv - 1);
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}

// [[Rcpp::export]]
List splitSE(
   NumericVector se, List ncat,
   int nroot, int nedge, int nleaf_unique,
   IntegerVector root, IntegerVector ulv, IntegerVector vlv,
   IntegerVector nclass, IntegerVector nclass_leaf
) {
   List lst_pi(nroot);
   List lst_tau(nedge);
   List lst_rho(nleaf_unique);

   double *se_ = se.begin();

   for (int r = 0; r < nroot; r ++) {
      NumericVector pi(nclass[root[r]] - 1);
      for (int k = 0; k < nclass[root[r]] - 1; k ++) {
         pi[k] = se_[k];
      }
      se_ += nclass[root[r]] - 1;
      lst_pi[r] = pi;
   }

   for (int d = 0; d < nedge; d ++) {
      int u = ulv[d]; int v = vlv[d];
      NumericMatrix tau(nclass[u] - 1, nclass[v]);
      for (int l = 0; l < nclass[v]; l ++) {
         for (int k = 0; k < nclass[u] - 1; k ++) {
            tau(k, l) = se_[k];
         }
         se_ += nclass[u] - 1;
      }
      lst_tau[d] = tau;
   }

   for (int v = 0; v < nleaf_unique; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector rho(nclass_leaf[v] * sum(ncatv - 1));
      double *rho_ = rho.begin();
      for (int k = 0; k < nclass_leaf[v]; k ++) {
         for (int m = 0; m < ncatv.length(); m ++) {
            for (int r = 0; r < ncatv[m] - 1; r ++) {
               rho_[r] = se_[r];
            }
            se_ += ncatv[m] - 1;
            rho_ += ncatv[m] - 1;
         }
      }
      lst_rho[v] = rho;
   }

   List res;
   res["pi"] = lst_pi;
   res["tau"] = lst_tau;
   res["rho"] = lst_rho;

   return res;
}
