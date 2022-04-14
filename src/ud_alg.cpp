#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
NumericVector pi_gnr(int nk, int nobs) {
   return rep(-log(nk), nk * nobs);
}

// [[Rcpp::export]]
NumericVector beta0_gnr(int nk, int ncov) {
   return rep(0.0, (nk - 1) * ncov);
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
NumericVector tau_gnr(int nl, int nk, int nobs) {
   NumericVector p(nk);
   NumericVector pv(nk * nl);
   double *ppv = pv.begin();

   for (int i = 0; i < nl; i ++) {
      p = runif(nk, 0, 1);
      for (int j = 0; j < nk; j ++) {
         ppv[j] = log(p[j] / sum(p));
      }
      ppv += nk;
   }

   return rep(pv, nobs);
}

// [[Rcpp::export]]
NumericVector beta1_gnr(int nk, int nl, int ncov) {
   return rep(0.0, (nk - 1) * nl * ncov);
}

// [[Rcpp::export]]
NumericVector beta2tau(
      NumericVector beta, NumericVector x,
      int nk, int nl, int np, int nobs
) {
   NumericVector pi(nobs * nk * nl);
   double *ppi = pi.begin();
   double *px = x.begin();
   for (int i = 0; i < nobs; i ++) {
      double *pb = beta.begin();
      for (int l = 0; l < nl; l ++) {
         double sl = 1;
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
      }
      px += np;
   }
   return pi;
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
   NumericVector pi;
   if (prob.isNull()) {
      pi = exp(pi_gnr(nk, 1));
   } else {
      pi = as<NumericVector>(prob);
   }

   IntegerVector cls = sample(nk, nobs, true, pi);
   return cls;
}

// [[Rcpp::export]]
IntegerVector cls_gnr(
   int nobs, int nk, int nl, IntegerVector vc,
   Nullable<NumericVector> prob = R_NilValue
) {
   NumericVector tau;
   if (prob.isNull()) {
      tau = exp(tau_gnr(nl, nk, 1));
   } else {
      tau = as<NumericVector>(prob);
   }

   NumericVector tau_l;
   IntegerVector cls(nobs);
   for (int i = 0; i < nobs; i ++) {
      tau_l = tau[seq_len(nk) + nk * (vc[i] - 1) - 1];
      cls[i] = sample(nk, 1, false, tau_l)[0];
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
      rho = exp(rho_gnr(nk, ncat));
   } else {
      rho = as<NumericVector>(prob);
   }

   int pos = 0;
   int nvar = ncat.length();
   NumericVector rho_m;
   IntegerVector y(nk * nvar);
   for (int i = 0; i < nobs; i ++) {
      pos = 0;
      int c = cls[i] - 1;
      for (int m = 0; m < nvar; m ++) {
         rho_m = rho[seq_len(ncat[m]) + pos + c * ncat[m] - 1];
         y[i * nvar + m] = sample(ncat[m], 1, false, rho_m)[0];
         pos += nk * ncat[m];
      }
   }

   return y;
}

// // [[Rcpp::export]]
// NumericVector y_gnr(
//    int nobs,
//    IntegerVector nvar, List ncat,
//    IntegerMatrix links, IntegerVector cstr_edge,
//    IntegerVector leafs, IntegerVector cstr_leaf,
//    int root,
//    IntegerVector nc,
//    Nullable<NumericVector> pi,
//    Nullable<List> tau,
//    Nullable<NumericVector> rho
// ) {
//    int nlv = nc.length();
//    if (pi.isNull()) NumericVector pi = pi_gnr(nc[root], 1);
//    if (tau.isNull()) {
//       List tau(links.nrow());
//       for (int d = 0; d < links.nrow(); d ++) {
//          NumericVector tau_d = tau_gnr(nc[links(d, 1)], nc[links(d, 0)], 1);
//          tau[d] = tau_d;
//       }
//    }
//    if (rho.isNull()) {
//       for (int v = 0; v < leafs.length(); v ++) {
//
//       }
//    }
//
//    List cls(nlv);
//    cls[root] = cls_gnr(nc[root], nobs, pi);
// }


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
   int nk, int nobs, int nvar, IntegerVector ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *rho = (double *)ptr_rho;
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

// nalpha = [#class] * nobs
// alpha = [#pclass] * nobs
// beta = [#pclass] * nobs
// ujbeta = [#pclass] * nobs
// lp = [#class * #pclass]
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
            val = exp(tau[k + l * nk] + ualpha[l] + ubeta[l] - jbeta[l]);
            joint[k + l * nk] += exp(val + beta[k] - ll[i]);
            svl += val;
         }
         alpha[k] = log(svl);
         post[k] += exp(alpha[k] + beta[k] - ll[i]);
      }
      tau += nk * nl;
      joint += nk * nl; post += nk;
      alpha += nk; ualpha += nl;
      beta += nk; ubeta += nl; jbeta += nl;
   }
}

void cumPosty(
   const double *denom, const double *numer, int *y,
   int nobs, int nvar, IntegerVector ncat, int nk,
   double *post, double *old_rho
) {
   for (int i = 0; i < nobs; i ++) {
      double *dnm = (double *)denom;
      double *nmr = (double *)numer;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            dnm[k] += exp(post[k]);
            if (y[m] > 0) nmr[y[m] - 1] += exp(post[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  nmr[r] += exp(post[k] + old_rho[r]);
               }
            }
            nmr += ncat[m];
            old_rho += ncat[m];
         }
         dnm += nk;
      }
      y += nvar;
   }
}

// y (nobs * sum(nvar))
// w (nobs)
// nvar (# measured lv)
// ncat (sum(nvar))
// links (nedge)
// constr (# measured lv)
// ncls (nlv)
// [[Rcpp::export]]
List treeFit(
   IntegerVector y, int nobs, IntegerVector nvar, List ncat,
   int nlv, int root, IntegerVector leaf, IntegerVector ulv, IntegerVector vlv,
   IntegerVector cstr_leaf, IntegerVector cstr_edge,
   IntegerVector nclass, IntegerVector ncleaf, IntegerVector nuclass, IntegerVector nvclass,
   int max_iter, double tol
) {
   int *py;
   int nleaf = leaf.length();
   int nedge = ulv.length();
   int nuleaf = ncleaf.length();
   int nuedge = nvclass.length();

   List lst_rho(nuleaf);
   List lst_tau(nuedge);
   List lst_rho_d(nuleaf);
   List lst_rho_n(nuleaf);
   List lst_ntau(nuedge);
   NumericVector pi = pi_gnr(nclass[root], nobs);
   std::vector<double*> ptr_rho(nuleaf);
   std::vector<double*> ptr_tau(nuedge);
   std::vector<double*> ptr_rho_d(nuleaf);
   std::vector<double*> ptr_rho_n(nuleaf);
   std::vector<double*> ptr_ntau(nuedge);

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

   for (int v = 0; v < nuleaf; v ++) {
      IntegerVector ncatv = ncat[v];
      NumericVector lrho = rho_gnr(ncleaf[v], ncatv);
      NumericVector rhod(ncleaf[v] * sum(ncatv));
      NumericVector rhon(ncleaf[v] * nvar[v]);
      ptr_rho[v] = lrho.begin();
      ptr_rho_d[v] = rhod.begin();
      ptr_rho_n[v] = rhon.begin();
      lst_rho[v] = lrho;
      lst_rho_d[v] = rhod;
      lst_rho_n[v] = rhon;
   }

   for (int d = 0; d < nuedge; d ++) {
      NumericVector ltau = tau_gnr(nuclass[d], nvclass[d], nobs);
      NumericVector ntau(nuclass[d] * nvclass[d] * nobs);
      ptr_tau[d] = ltau.begin();
      ptr_ntau[d] = ntau.begin();
      lst_tau[d] = ltau;
      lst_ntau[d] = ntau;
   }

   for (int d = 0; d < nedge; d ++) {
      int lk = cstr_edge[d];
      NumericMatrix joint(nuclass[lk] * nvclass[lk], nobs);
      ptr_joint[d] = joint.begin();
      lst_joint[d] = joint;
      NumericVector jbeta(nvclass[lk] * nobs);
      ptr_j[d] = jbeta.begin();
      lst_j[d] = jbeta;
   }

   for (int v = 0; v < nlv; v ++) {
      NumericVector post(nclass[v]);
      NumericVector alpha(nclass[v] * nobs);
      NumericVector beta(nclass[v] * nobs);
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

      // (expectation-step)
      // initiate beta
      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         msrPrd(py, ptr_rho[cstr_leaf[v]], ptr_b[leaf[v]], ncleaf[cstr_leaf[v]],
                nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]]);
         py += nvar[cstr_leaf[v]] * nobs;
      }

      // upward recursion
      for (int d = 0; d < nedge; d ++) {
         int u = ulv[d];
         int v = vlv[d];
         upRec(ptr_b[v], ptr_j[d], ptr_b[u], ptr_tau[cstr_edge[d]],
               nobs, nclass[u], nclass[v]);
      }

      // initiate alpha
      lst_a[root] = pi;
      double *beta1 = ptr_b[root];
      double *alpha1 = ptr_a[root];
      double *post1 = ptr_post[root];
      for (int i = 0; i < nobs; i ++) {
         double lik = 0;
         for (int k = 0; k < nclass[root]; k ++) {
            post1[k] = alpha1[k] + beta1[k];
            lik += exp(post1[k]);
         }
         ll[i] = log(lik);
         currll += ll[i];

         for (int k = 0; k < nclass[root]; k ++) {
            post1[k] -= ll[i];
         }

         alpha1 += nclass[root];
         beta1  += nclass[root];
         post1  += nclass[root];
      }

      // Downward recursion
      for (int d = nedge - 1; d == 0; d --) {
         int u = ulv[d];
         int v = vlv[d];
         dnRec(ptr_a[u], ptr_a[v], ptr_b[u], ptr_b[v], ptr_j[d],
               nobs, nclass[u], nclass[v], ptr_tau[cstr_edge[d]],
               ptr_post[u], ptr_joint[d], ll);
      }

      // (maximization-step)
      // pi updates
      NumericVector new_pi(nclass[root]);
      double *post = ptr_post[root];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < nclass[root]; k ++) {
            new_pi[k] += post[k];
         }
         post += nclass[root];
      }
      pi = rep(log(new_pi / sum(new_pi)), nobs);

      // tau updates
      for (int d = 0; d < nedge; d ++) {
         NumericVector ntau = lst_ntau[cstr_edge[d]];
         NumericMatrix joint = lst_joint[d];
         NumericVector jsum = rowSums(joint);
         ntau += jsum;
      }
      for (int d = 0; d < nuedge; d ++) {
         double *ptau = ptr_ntau[d];
         for (int l = 0; l < nvclass[d]; l ++) {
            double sl = 0;
            for (int k = 0; k < nuclass[d]; k ++) {
               sl += exp(ptau[k]);
            }
            sl = log(sl);
            for (int k = 0; k < nuclass[d]; k ++) {
               ptau[k] -= sl;
            }
            ptau += nuclass[d];
         }
         NumericVector ntau = lst_ntau[d];
         lst_tau[d] = rep(ntau, nobs);
      }

      py = y.begin();
      for (int v = 0; v < nleaf; v ++) {
         cumPosty(ptr_rho_d[cstr_leaf[v]], ptr_rho_n[cstr_leaf[v]],
                  py, nobs, nvar[cstr_leaf[v]], ncat[cstr_leaf[v]],
                  ncleaf[v], ptr_post[v], ptr_rho[cstr_leaf[v]]);
         py += nobs * nvar[v];
      }
      for (int v = 0; v < nuleaf; v ++) {
         double *rho = ptr_rho[v];
         double *numer = ptr_rho_n[v];
         double *denom = ptr_rho_d[v];
         IntegerVector ncatv = ncat[v];
         for (int m = 0; m < nvar[v]; m ++) {
            for (int k = 0; k < ncleaf[v]; k ++) {
               for (int r = 0; r < ncatv[m]; r ++) {
                  rho[r] = log(numer[r] / denom[k]);
               }
               rho  += ncatv[m];
               numer += ncatv[m];
            }
            denom += ncleaf[v];
         }
      }

      dll = (currll - lastll) / lastll;
      Rcout << iter << ") ll: " << currll << " / diff: " << dll << "\n" << std::endl;
   }

   // computes posterior probs
   double loglik = currll;
   List posterior(nlv);
   for (int v = 0; v < nlv; v ++) {
      NumericMatrix post(nobs, nclass[v]);
      double *alpha = ptr_a[v];
      double *beta  = ptr_b[v];
      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < nclass[v]; k ++) {
            post(i, k) = exp(alpha[k] + beta[k] - ll[i]);
         }
         alpha += nclass[v];
         beta  += nclass[v];
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
NumericMatrix test1(
      NumericMatrix a, NumericMatrix b
) {
   a += b;
   return a;
}
