#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
int sample1(int n, NumericVector prob) {
   double r = R::runif(0, 1);
   double cp = 0;
   for (int i = 0; i < prob.length() - 1; i ++) {
      cp += prob[i];
      if (r < cp) return i;
   }
   return n - 1;
}

// [[Rcpp::export]]
NumericMatrix pi_gnr(int nk, int nobs) {
   NumericMatrix pi(nk, nobs);
   pi.fill(-log(nk));
   return pi;
}

// [[Rcpp::export]]
NumericMatrix beta0_gnr(int nk, int p) {
   NumericMatrix beta(nk - 1, p);
   beta.fill(0);
   return beta;
}

// [[Rcpp::export]]
NumericVector beta2pi(
   NumericVector beta, NumericVector x,
   int nk, int np, int nobs
) {
   NumericVector pi(nobs * nk);
   double *ppi = pi.begin();
   double *px = x.begin();
   for (int i = 0; i < nobs; i ++) {
      double sl = 1;
      double *pb = beta.begin();
      for (int k = 0; k < nk - 1; k ++) {
         for (int p = 0; p < np; p ++) {
            ppi[k] += pb[p] * px[p];
         }
         sl += exp(ppi[k]);
         pb += np;
      }

      for (int k = 0; k < nk; k ++) {
         ppi[k] -= log(sl);
      }
      ppi += nk;
      px  += np;
   }
   return pi;
}

// [[Rcpp::export]]
NumericMatrix tau_gnr(int nk, int nl, int nobs) {
   NumericMatrix tau(nk * nl, nobs);
   NumericMatrix ntau(nk, nl);
   double *ptau = tau.begin();

   for (int l = 0; l < nl; l ++) {
      NumericVector p = runif(nk, 0, 1);
      ntau.column(l) = log(p) - log(sum(p));
   }
   for (int i = 0; i < nobs; i ++) {
      double *pntau = ntau.begin();
      for (int l = 0; l < nl; l ++) {
         for (int k = 0; k < nk; k ++)
            ptau[k] = pntau[k];
         pntau += nk;
         ptau  += nk;
      }
   }

   return tau;
}

// [[Rcpp::export]]
NumericVector beta1_gnr(int nk, int nl, int p) {
   return rnorm((nk - 1) * nl * p, 0, 0.1);
}

// [[Rcpp::export]]
NumericMatrix beta2tau(
   NumericVector beta, NumericMatrix x,
   int nk, int nl, int np, int nobs
) {
   NumericMatrix tau(nk * nl, nobs);
   double *ptau = tau.begin();
   double *px = x.begin();
   for (int i = 0; i < nobs; i ++) {
      double *pbeta = beta.begin();
      for (int l = 0; l < nl; l ++) {
         double sl = 1;
         for (int k = 0; k < nk - 1; k ++) {
            for (int p = 0; p < np; p ++) {
               ptau[k] += pbeta[p] * px[p];
            }
            sl += exp(ptau[k]);
            pbeta += np;
         }
         for (int k = 0; k < nk; k ++) {
            ptau[k] -= log(sl);
         }
         ptau += nk;
      }
      px += np;
   }
   return tau;
}

// [[Rcpp::export]]
NumericVector rho_gnr(int nk, IntegerVector ncat) {
   NumericVector pv(nk * sum(ncat));
   double *ppv = pv.begin();

   for (int m = 0; m < ncat.length(); m ++) {
      for (int i = 0; i < nk; i ++) {
         NumericVector p = runif(ncat[m], 0, 1);
         for (int j = 0; j < ncat[m]; j ++) {
            ppv[j] = log(p[j] / sum(p));
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
   NumericMatrix pi;
   if (prob.isNull()) {
      pi = pi_gnr(nk, nobs);
   } else {
      pi = as<NumericMatrix>(prob);
   }

   IntegerVector cls(nobs);
   for (int i = 0; i < nobs; i ++) {
      cls[i] = sample1(nk, exp(pi.column(i)));
   }

   return cls;
}

// [[Rcpp::export]]
IntegerVector cls_gnr(
   int nobs, int nk, int nl, IntegerVector v,
   Nullable<NumericVector> prob = R_NilValue
) {
   NumericMatrix tau;
   if (prob.isNull()) {
      tau = tau_gnr(nl, nk, nobs);
   } else {
      tau = as<NumericMatrix>(prob);
   }

   IntegerVector cls(nobs);
   double *ptau = tau.begin();
   for (int i = 0; i < nobs; i ++) {
      NumericMatrix tau_l(nk, nl, ptau);
      cls[i] = sample1(nk, exp(tau_l.column(v[i])));
      ptau += nk * nl;
   }

   return cls;
}

// [[Rcpp::export]]
IntegerVector y_gnr(
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
   NumericVector rho_m;
   IntegerVector y(nobs * nvar);
   for (int i = 0; i < nobs; i ++) {
      double *pos = rho.begin();
      for (int m = 0; m < nvar; m ++) {
         NumericMatrix rho_m(ncat[m], nk, pos);;
         y[i * nvar + m] = sample1(ncat[m], exp(rho_m.column(cls[i]))) + 1;
         pos += nk * ncat[m];
      }
   }

   return y;
}

// Upward-Recursion
void upInit(
   int *y, const double *ptr_rho, double *beta,
   int nk, int nobs, int nvar, IntegerVector ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *rho = (double *) ptr_rho;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0)
               beta[k] += rho[y[m] - 1];
            rho += ncat[m];
         }
      }
      y += nvar;
      beta += nk;
   }
}

void upRec(
   double *beta, double *jbeta, double *lbeta,
   const double *tau, int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(tau[k] + lbeta[k]);

         tau += nk;
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      lbeta += nk;
      jbeta += nl;
      beta  += nl;
   }
}

// Downward-Recursion
void dnInit(
   double *alpha, double *beta, double *pi,
   double *post, double *ll, int nobs, int nclass
) {
   for (int i = 0; i < nobs; i ++) {
      double lik = 0;
      for (int k = 0; k < nclass; k ++) {
         alpha[k] = pi[k];
         post[k] = alpha[k] + beta[k];
         lik += exp(post[k]);
      }
      ll[i] += log(lik);

      for (int k = 0; k < nclass; k ++) {
         post[k] -= ll[i];
      }

      pi    += nclass;
      alpha += nclass;
      beta  += nclass;
      post  += nclass;
   }
}

void dnRec(
   double *alpha, double *ualpha,
   double *beta, double *ubeta, double *jbeta,
   int nobs, int nk, int nl, double *tau,
   double *post, double *joint, NumericVector ll
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         double svl = 0;
         for (int l = 0; l < nl; l ++) {
            val = tau[k + l * nk] + ualpha[l] + ubeta[l] - jbeta[l];
            joint[k + l * nk] = val + beta[k] - ll[i];
            svl += exp(val);
         }
         alpha[k] = log(svl);
         post[k] = alpha[k] + beta[k] - ll[i];
      }
      tau += nk * nl;
      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      beta += nk; ubeta += nl; jbeta += nl;
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
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nclass; k ++) {
         pi[k] = log(npi[k]) - log(sum(npi));
      }
      pi += nclass;
   }
}

void cumTau(
   double *joint, const double *ptr_ntau,
   int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *ntau = (double *) ptr_ntau;
      for (int l = 0; l < nl; l ++) {
         for (int k = 0; k < nk; k ++) {
            ntau[k] += exp(joint[k]);
         }
         joint += nk;
         ntau  += nk;
      }
   }
}

void updateTau(
   double *tau, const double *ptr_ntau,
   int nobs, int nk, int nl
) {
   for (int i = 0; i < nobs; i ++) {
      double *ntau = (double *) ptr_ntau;
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
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            dnm[k] += exp(post[k]);
            if (y[m] > 0) nmr[y[m] - 1] += exp(post[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  nmr[r] += exp(post[k] + rho[r]);
               }
            }
            nmr += ncat[m];
            rho += ncat[m];
         }
         dnm += nk;
      }
      post += nk;
      y += nvar;
   }
}

void updateRho(
   double *rho, double *numer, double *denom,
   int nobs, int nclass, int nvar, IntegerVector ncat
) {
   for (int m = 0; m < nvar; m ++) {
      for (int k = 0; k < nclass; k ++) {
         for (int r = 0; r < ncat[m]; r ++) {
            rho[r] = log(numer[r]) - log(denom[k]);
         }
         rho   += ncat[m];
         numer += ncat[m];
      }
      denom += nclass;
   }
}

// [[Rcpp::export]]
List treeFit(
      IntegerVector y, int nobs, IntegerVector nvar, List ncat,
      int nlv, IntegerVector root, IntegerVector leaf,
      IntegerVector ulv, IntegerVector vlv,
      IntegerVector cstr_leaf, IntegerVector cstr_edge,
      IntegerVector nclass, IntegerVector nclass_leaf,
      IntegerVector nclass_u, IntegerVector nclass_v,
      LogicalVector init, List init_param,
      int max_iter, double tol
) {
   int *py;
   int nroot = root.length();
   int nleaf = leaf.length();
   int nedge = ulv.length();
   int nleaf_unique = nclass_leaf.length();
   int nedge_unique = nclass_v.length();

   List lst_pi(nroot);
   List lst_tau(nedge_unique);
   List lst_rho(nleaf_unique);
   List lst_ntau(nedge_unique);
   List lst_nrho_d(nleaf_unique);
   List lst_nrho_n(nleaf_unique);
   std::vector<double*> ptr_pi(nroot);
   std::vector<double*> ptr_tau(nedge_unique);
   std::vector<double*> ptr_rho(nleaf_unique);
   std::vector<double*> ptr_ntau(nedge_unique);
   std::vector<double*> ptr_nrho_d(nleaf_unique);
   std::vector<double*> ptr_nrho_n(nleaf_unique);

   List lst_a(nlv);
   List lst_b(nlv);
   List lst_j(nedge);
   std::vector<double*> ptr_a(nlv);
   std::vector<double*> ptr_b(nlv);
   std::vector<double*> ptr_j(nedge);

   List lst_post(nlv);
   List lst_joint(nedge);
   NumericVector ll(nobs);
   std::vector<double*> ptr_post(nlv);
   std::vector<double*> ptr_joint(nedge);

   if (init[0]) {
      List piList = init_param["pi"];
      for (int r = 0; r < nroot; r ++) {
         NumericMatrix pi_init = piList[r];
         NumericMatrix pi = clone(pi_init);
         ptr_pi[r] = pi.begin();
         lst_pi[r] = pi;
      }
   } else {
      for (int r = 0; r < nroot; r ++) {
         NumericMatrix pi = pi_gnr(nclass[root[r]], nobs);
         ptr_pi[r] = pi.begin();
         lst_pi[r] = pi;
      }
   }

   if (init[1]) {
      List tauList = init_param["tau"];
      for (int d = 0; d < nedge_unique; d ++) {
         NumericMatrix tau_init = tauList[d];
         NumericMatrix tau = clone(tau_init);
         NumericMatrix ntau(nclass_u[d], nclass_v[d]);
         ptr_tau[d] = tau.begin();
         ptr_ntau[d] = ntau.begin();
         lst_tau[d] = tau;
         lst_ntau[d] = ntau;
      }
   } else {
      for (int d = 0; d < nedge_unique; d ++) {
         NumericMatrix tau = tau_gnr(nclass_u[d], nclass_v[d], nobs);
         NumericMatrix ntau(nclass_u[d], nclass_v[d]);
         ptr_tau[d] = tau.begin();
         ptr_ntau[d] = ntau.begin();
         lst_tau[d] = tau;
         lst_ntau[d] = ntau;
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

   for (int d = 0; d < nedge; d ++) {
      int lk = cstr_edge[d];
      NumericMatrix joint(nclass_u[lk] * nclass_v[lk], nobs);
      ptr_joint[d] = joint.begin();
      lst_joint[d] = joint;
      NumericMatrix jbeta(nclass_v[lk], nobs);
      ptr_j[d] = jbeta.begin();
      lst_j[d] = jbeta;
   }

   for (int v = 0; v < nlv; v ++) {
      NumericMatrix post(nclass[v], nobs);
      NumericMatrix alpha(nclass[v], nobs);
      NumericMatrix beta(nclass[v], nobs);
      ptr_post[v] = post.begin();
      ptr_a[v] = alpha.begin();
      ptr_b[v] = beta.begin();
      lst_post[v] = post;
      lst_a[v] = alpha;
      lst_b[v] = beta;
   }


   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;

   while ( (iter < max_iter) && (dll > tol) ) {
      iter ++;
      lastll = currll;

      // beta, cleaning
      for (int v = 0; v < nlv; v ++) {
         NumericMatrix beta = lst_b[v];
         beta.fill(0);
      }

      // (expectation-step)
      // initiate beta
      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         upInit(py, ptr_rho[cstr_leaf[v]], ptr_b[leaf[v]], nclass[leaf[v]],
                nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
         py += nobs * nvar[cstr_leaf[v]];
      }

      // upward recursion
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d];
         int v = vlv[d];
         upRec(ptr_b[v], ptr_j[d], ptr_b[u], ptr_tau[cstr_edge[d]],
               nobs, nclass[u], nclass[v]);
      }

      // initiate alpha
      ll.fill(0);
      for (int r = 0; r < nroot; r ++) {
         dnInit(ptr_a[root[r]], ptr_b[root[r]], ptr_pi[r], ptr_post[root[r]],
                ll.begin(), nobs, nclass[root[r]]);
      }

      // Downward recursion
      for (int d = nedge - 1; d > -1; d --) {
         int u = ulv[d];
         int v = vlv[d];
         dnRec(ptr_a[u], ptr_a[v], ptr_b[u], ptr_b[v], ptr_j[d],
               nobs, nclass[u], nclass[v], ptr_tau[cstr_edge[d]],
               ptr_post[u], ptr_joint[d], ll);
      }

      // (maximization-step)
      // pi updates
      for (int r = 0; r < nroot; r ++) {
         updatePi(ptr_pi[r], ptr_post[root[r]],
                  nobs, nclass[root[r]]);
      }

      // tau updates
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d]; int v = vlv[d];
         cumTau(ptr_joint[d], ptr_ntau[cstr_edge[d]],
                nobs, nclass[u], nclass[v]);
      }
      for (int d = 0; d < nedge_unique; d ++) {
         int nk = nclass_u[d]; int nl = nclass_v[d];
         updateTau(ptr_tau[d], ptr_ntau[d], nobs, nk, nl);
      }

      // rho updates
      py = y.begin();
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

      currll = sum(ll);
      if (lastll == R_NegInf) dll = R_PosInf;
      else dll = currll - lastll;
      dll = 1;
      Rcout << iter << ") ll: " << currll << " / diff: " << dll << std::endl;
   }

   // computes posterior probs
   double loglik = currll;
   List posterior(nlv);
   for (int v = 0; v < nlv; v ++)
      posterior[v] = lst_post[v];

   List par;
   par["pi"] = lst_pi;
   par["tau"] = lst_tau;
   par["rho"] = lst_rho;

   List ret;
   ret["param"] = par;
   ret["posterior"] = posterior;
   ret["loglik"] = loglik;

   return ret;
}


// beta is already included prevalence
void upRecS(
   double *nbeta, double *jbeta, double *beta,
   double *prv, double *lp,
   int nlv, int nl, int *nk
) {
   double nu = 0;
   for (int j = 0; j < nlv; j ++) {
      for (int l = 0; l < nl; l ++) {
         double mll = 0;
         for (int k = 0; k < nk[j]; k ++)
            mll += exp(lp[k] + beta[k] - prv[k]);

         lp += nk[j];
         jbeta[l] = log(mll);
         nbeta[l] += log(mll);
      }

      jbeta += nl;
      beta += nk[j];
   }

   for (int l = 0; l < nlv; l ++)
      nu += exp(nbeta[l]);
   for (int l = 0; l < nlv; l ++)
      nbeta[l] -= nu;
}


void dnRecS(
   double *pst, double *jpst, double *upst,
   double *beta, double *jbeta,
   double *prv, double *lp,
   int nl, int nk
) {
   for (int k = 0; k < nk; k ++) {
      double fwd = beta[k] - prv[k];
      double sbwd = 0;
      for (int l = 0; l < nl; l ++) {
         double bwd = lp[k + l * nk] + upst[l] - jbeta[l];
         jpst[k + l * nk] = fwd + bwd;
         sbwd += exp(bwd);
      }
      pst[k] = fwd + log(sbwd);
   }
}

// [[Rcpp::export]]
NumericMatrix test1(
      NumericMatrix a, NumericMatrix b
) {
   a += b;
   return a;
}
