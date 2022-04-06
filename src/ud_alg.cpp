#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector pi_gnr(int k, int n) {
   NumericVector p(k);
   p = runif(k, 0, 1);
   return rep(log(p / sum(p)), n);
}

// [[Rcpp::export]]
NumericVector tau_gnr(int l, int k, int n) {
   NumericVector p(k);
   NumericVector pv(k * l);
   double *ppv = pv.begin();

   for (int i = 0; i < l; i ++) {
      p = runif(k, 0, 1);
      for (int j = 0; j < k; j ++) {
         ppv[j] = log(p[j] / sum(p));
      }
      ppv += k;
   }

   return rep(pv, n);
}


// [[Rcpp::export]]
NumericVector rho_gnr(int k, IntegerVector ncat) {
   NumericVector pv(k * sum(ncat));
   double *ppv = pv.begin();

   for (int m = 0; m < ncat.length(); m ++) {
      for (int i = 0; i < k; i ++) {
         NumericVector p = runif(ncat[m], 0, 1);
         for (int j = 0; j < ncat[m]; j ++) {
            ppv[j] = log(p[j] / sum(p));
         }
         ppv += ncat[m];
      }
   }
   return pv;
}


// measured lv : beta[nk * obs] / alpha[nk * obs]
// middle lv : beta[nl * obs] / jbeta[nk * nl * obs] / alpha[nk * obs]
// summit lv : beta[nl * obs] / alpha[nl * obs]

double logAdd(double lx, double ly) {
   if (lx == R_NegInf) return ly;
   if (ly == R_NegInf) return lx;

   double lxy = 0;
   if (lx < ly)
      lxy = ly + log(1 + exp(lx - ly));
   else
      lxy = lx + log(1 + exp(ly - lx));

   return lxy;
}


void msrPrd(
   int *y, const double *ptr_rho, double *beta,
   int nk, int nobs, int nvar, int *ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *ptr_rhoi = (double *)ptr_rho;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0)
               beta[k] += ptr_rhoi[y[m] - 1];

            ptr_rhoi += ncat[m];
         }
      }
      y += nvar;
      beta += nk;
   }
}

void argMpr(
   int *y, double *pst,
   const double *ptr_rho, const double *npr,
   int nk, int nobs, int nvar, int *ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *ptr_rhoi = (double *) ptr_rho;
      double *npri = (double *) npr;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0) npri[y[m] - 1] += exp(pst[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  npri[r] += exp(pst[k] + ptr_rho[r]);
               }
            }
            y += nvar;
            ptr_rhoi += ncat[m];
            npri += ncat[m];
         }
      }
      pst += nk;
   }
}


// beta = [#class] * nobs
// jbeta = [#class * nlv] * nobs
// lst_b = [#cclass per nlv] * nobs
// lp = [#class * #pclass]
void upRec(
   double *beta, double *jbeta, double *lbeta,
   double *tau, int nobs, int nl, int nk
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
      jbeta += nl;
      beta  += nl;
   }
}

void upRec2(
      double *beta, double *jbeta, double *lst_b,
      double *ltau, int nobs, int nl, int nk
) {
   for (int j = 0; j < nobs; j ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(ltau[k] + lst_b[k]);

         ltau += nk;
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
      lst_b += nk;
   }
}

// nalpha = [#class] * nobs
// alpha = [#pclass] * nobs
// beta = [#pclass] * nobs
// ujbeta = [#pclass] * nobs
// lp = [#class * #pclass]
void dnRec(
   double *alpha, double *ualpha,
   double *beta, double *ubeta, double *jbeta,
   int nobs, int nl, int nk, double *tau,
   double *post, double *joint, NumericVector ll
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         double svl = 0;
         for (int l = 0; l < nl; l ++) {
            val = exp(tau[k + l * nk] + ualpha[l] + ubeta[l] - jbeta[l]);
            joint[k + l * nk] += exp(val + beta[k] - ll[i]);
            svl += val;
         }
         alpha[k] = log(svl);
         post[k] += exp(alpha[k] + beta[k] - ll[i]);
      }
      tau += nk * nl;
      alpha += nk; ualpha += nl;
      beta += nk; ubeta += nl; jbeta += nl;
   }
}

// y (nobs * sum(nvar))
// w (nobs)
// nvar (# measured lv)
// ncat (sum(nvar))
// links (nlink)
// constr (# measured lv)
// ncls (nlv)
List treeFit(
   IntegerVector y, NumericVector logw,
   IntegerVector nvar, IntegerVector ncat,
   IntegerVector leaf, IntegerVector cstr_lf,
   IntegerMatrix links, IntegerVector cstr_lk,
   IntegerVector ncls,
   IntegerVector lnc, IntegerVector unc,
   int max_iter, double tol
) {
   int nobs = logw.length();
   int nlv  = ncls.length();
   int nleaf = leaf.length();
   int nlink = links.nrow();
   int nuleaf = max(cstr_lf);
   int nulink = max(cstr_lk);

   List lst_rho(nuleaf);
   List lst_tau(nulink);
   NumericVector lw = rep_each(logw, ncls[nlv - 1]);
   NumericVector pi = pi_gnr(ncls[nlv - 1], nobs) + lw;
   std::vector<double*> ptr_rho(nleaf);
   std::vector<double*> ptr_tau(nlink);

   List lst_a(nlv);
   List lst_b(nlv);
   List lst_j(nlink);
   std::vector<double*> ptr_a(nlv);
   std::vector<double*> ptr_b(nlv);
   std::vector<double*> ptr_j(nlink);

   List lst_post(nlv);
   List lst_joint(nlink);
   NumericVector ll(nobs);
   std::vector<double*> ptr_post(nlv);
   std::vector<double*> ptr_joint(nlink);

   for (int v = 0; v < nuleaf; v ++) {
      NumericVector lrho = rho_gnr(ncls[v], ncat[v]);
      ptr_rho[v] = lrho.begin();
      lst_rho[v] = lrho;
   }

   for (int v = 0; v < nleaf; v ++) {
      NumericVector post(ncls[v]);
      ptr_post[v] = post.begin();
      lst_post[v] = post;
   }

   for (int d = 0; d < nulink; d ++) {
      NumericVector ltau = tau_gnr(lnc[d], unc[d], nobs);
      ptr_tau[d] = ltau.begin();
      lst_tau[d] = ltau;
   }

   for (int d = 0; d < nlink; d ++) {
      NumericVector joint(ncls[links(d, 0)] * ncls[links(d, 1)]);
      ptr_joint[d] = joint.begin();
      lst_joint[d] = joint;
      NumericVector jbeta(nobs * ncls[links(d, 1)]);
      ptr_j[d] = jbeta.begin();
      lst_j[d] = jbeta;
   }

   for (int v = 0; v < nlv; v ++) {
      NumericVector alpha(nobs * ncls[v]);
      NumericVector beta(nobs * ncls[v]);
      ptr_a[v] = alpha.begin();
      ptr_b[v] = beta.begin();
      lst_a[v] = alpha;
      lst_b[v] = beta;
   }

   int iter = 0;
   double currll = R_NegInf;
   double lastll = R_NegInf;
   double dll = R_PosInf;
   while (iter < max_iter || (dll > tol)) {
      iter ++;
      lastll = currll;

      // (expectation-step)
      // initiate beta
      int *py = y.begin();
      int *pnc = ncat.begin();
      for (int v = 0; v < nleaf; v ++) {
         msrPrd(py, ptr_rho[cstr_lf[v]], ptr_b[v],
                ncls[v], nobs, nvar[v], pnc);
         py += nvar[v] * nobs;
         pnc += nvar[v];
      }

      // upward recursion
      for (int d = 0; d < nlink; d ++) {
         int u = links(d, 0);
         int v = links(d, 1);
         upRec(ptr_b[v], ptr_j[d], ptr_b[u], ptr_tau[cstr_lk[d]],
               nobs, ncls[u], ncls[v]);
      }

      // initiate alpha
      lst_a[nlv - 1] = pi;
      double *beta1 = ptr_b[nlv - 1];
      double *alpha1 = ptr_a[nlv - 1];
      double *post1 = ptr_post[nlv - 1];
      for (int i = 0; i < nobs; i ++) {
         double lik = 0;
         for (int k = 0; k < ncls[nlv - 1]; k ++) {
            post1[k] = alpha1[k] + beta1[k];
            lik += exp(post1[k]);
         }
         ll[i] = log(lik);
         currll += ll[i];

         for (int k = 0; k < ncls[nlv - 1]; k ++) {
            post1[k] -= ll[i];
         }

         alpha1 += ncls[nlv - 1];
         beta1  += ncls[nlv - 1];
         post1  += ncls[nlv - 1];
      }

      // Downward recursion
      for (int d = nlink - 1; d == 0; d --) {
         int u = links(d, 0);
         int v = links(d, 1);
         int nl = ncls[v];  int nk = ncls[u];
         dnRec(ptr_a[u], ptr_a[v], ptr_b[u], ptr_b[v], ptr_j[d],
               nobs, ncls[v], ncls[u], ptr_tau[cstr_lk[d]],
               ptr_post[cstr_lk[d]], ptr_joint[cstr_lk[d]], ll);
      }

      // (maximization-step)
      // pi updates
      NumericVector new_pi(ncls[nlv - 1]);
      double *post = ptr_post[nlv - 1];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < ncls[nlv - 1]; k ++) {
            new_pi[k] += post[k];
         }
         post += ncls[nlv - 1];
      }
      pi = rep(new_pi / sum(new_pi), nobs);

      // tau updates
      std::vector<double*> denom(nulink);
      std::vector<double*> numer(nulink);
      for (int d = 0; d < nulink; d++) {
         NumericVector den(lnc[d] * unc[d]);
         NumericVector num(lnc[d]);
         denom[d] = den.begin();
         numer[d] = num.begin();
      }
      for (int d = 0; d < nlink; d ++) {
         double *post  = ptr_post[d];
         double *joint = ptr_joint[d];
         for (int i = 0; i < nobs; i ++) {
            for (int l = 0; l < ncls[links(d, 1)]; l ++) {
               for (int k = 0; k < ncls[links(d, 0)]; k ++) {

               }
            }
         }
      }

      List nlst_rho(nleaf);
      // for (int v = 0; v < nleaf; v ++) {
      //    NumericVector opr = lst_rho[v];
      //    NumericVector npr(opr.length());
      //    int *py  = y.begin();
      //    double *rho = ptr_rho[cstr_lf[v]];
      //    double *pnpr = npr.begin();
      //    double *marg = ptr_marg[v];
      //    for (int m = 0; m < nvar[v]; m ++) {
      //       for (int k = 0; k < ncls[v]; k ++) {
      //          for (int i = 0; i < nobs; i ++) {
      //             if (py[i] > 0) ptr_rho[py[i]] += exp(marg[v]);
      //             else {
      //                for (int r = 0; r < ncat[m]; r ++) {
      //                   pnpr[r] += exp(marg[v] + ptr_rho[r]);
      //                }
      //             }
      //          }
      //          py += nobs;
      //       }
      //       pnpr += ncat[m];
      //       ptr_rho += ncat[m];
      //    }
      // }
   }

   // computes posterior probs
   double loglik = currll;
   List posterior(nlv);
   for (int v = 0; v < nlv; v ++) {
      NumericMatrix post(nobs, ncls[v]);
      double *alpha = ptr_a[v];
      double *beta  = ptr_b[v];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < ncls[v]; k ++) {
            post(i, k) = exp(alpha[k] + beta[k] - ll[i]);
         }
         alpha += ncls[v];
         beta  += ncls[v];
      }
      posterior[v] = post;
   }

   List ret;
   ret["pi"] = pi;
   ret["tau"] = lst_tau;
   ret["rho"] = lst_rho;
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
List test1(
      int n, IntegerVector k
) {
   List l(n);
   std::vector<int*> v(0);

   for (int i = 0; i < n; i ++) {
      IntegerVector li(k[i] * k[i]);
      v.push_back(li.begin());
      l[i] = li;
   }
   for (int i = 0; i < n; i ++) {
      int *vi = v[i];
      for (int j = 0; j < k[i]; j ++) {
         for (int m = 0; m < k[i]; m ++) {
            vi[m] = i * j * m;
         }
         vi += k[i];
      }
   }
   return l;
}
