#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// class catlv{
// public:
//    int nk;
//
//    int depth;
//    bool isleaf;
//    bool isroot;
//
//    int parent;
//    int *children;
//
//    double *beta;
//    double *alpha;
//
//    int *y;
//
//
//    catlv(List _info) {
//       nk = _info["nclass"];
//       depth = _info["depth"];
//       isleaf = _info["isleaf"];
//       isroot = _info["isroot"];
//       parent = _info["parent"];
//       children = _info["children"];
//       y = _info["y"];
//    }
//
//    // void msrPrd();
//    // void upRec();
//    // void dnRec();
//    // void udPrd();
//    // void upRecS();
//    // void dnRecS();
//    // void udPrdS();
// };

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
   int *y, const double *prho, double *beta,
   int nk, int nobs, int nvar, int *ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *prhoi = (double *)prho;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0)
               beta[k] += prhoi[y[m] - 1];

            prhoi += ncat[m];
         }
      }
      y += nvar;
      beta += nk;
   }
}

void argMpr(
   int *y, double *pst,
   const double *prho, const double *npr,
   int nk, int nobs, int nvar, int *ncat
) {
   for (int i = 0; i < nobs; i ++) {
      double *prhoi = (double *) prho;
      double *npri = (double *) npr;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0) npri[y[m] - 1] += exp(pst[k]);
            else {
               for (int r = 0; r < ncat[m]; r ++) {
                  npri[r] += exp(pst[k] + prho[r]);
               }
            }
            y += nvar;
            prhoi += ncat[m];
            npri += ncat[m];
         }
      }
      pst += nk;
   }
}


// beta = [#class] * nobs
// jbeta = [#class * nlvar] * nobs
// lb = [#cclass per nlvar] * nobs
// lp = [#class * #pclass]
void upRec(
   double *beta, double *jbeta, double *lb,
   double *prt, int nobs, int nl, int nk
) {
   for (int i = 0; i < nobs; i ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(prt[k] + lb[k]);

         prt += nk;
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
      beta  += nl;
   }
}

void upRec2(
      double *beta, double *jbeta, double *lb,
      double *lprt, int nobs, int nl, int nk
) {
   for (int j = 0; j < nobs; j ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(lprt[k] + lb[k]);

         lprt += nk;
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
      lb += nk;
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
   int nobs, int nl, int nk, double *prt,
   double *mpt, double *umpt, double *jpt
) {
   for (int i = 0; i < nobs; i ++) {
      for (int k = 0; k < nk; k ++) {
         double val = 0;
         double svl = 0;
         for (int l = 0; l < nl; l ++) {
            val = exp(prt[k + l * nk] + umpt[l] - jbeta[l]);
            jpt[k + l * nk] = val + beta[k];
            svl += val;
         }
         alpha[k] = log(svl);
         mpt[k] = alpha[k] + beta[k];
      }
      prt += nk * nl;
      jpt += nk * nl; mpt += nk;
      alpha += nk; ualpha += nl;
      beta += nk; ubeta += nl; jbeta += nl;
   }
}

void edgePrb(
      double *joint, double *ualpha,
      double *beta, double *jbeta, double *ubeta,
      double *lprt, int nl, int nk
) {
   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk; k ++)
         joint[l] = lprt[k] + ualpha[l] +
            beta[k] + ubeta[l] - jbeta[l];

      lprt += nk;
      joint += nk;
   }
}

// y (nobs * sum(nvar))
// w (nobs)
// nvar (# measured lv)
// ncat (sum(nvar))
// edges (nedge)
// constr (# measured lv)
// ncls (nlvar)
List treeFit(
   IntegerVector y, NumericVector w,
   IntegerVector leaf, IntegerVector constr,
   IntegerVector nvar, IntegerVector ncat,
   IntegerMatrix edges, IntegerVector ncls,
   int max_iter, double tol
) {
   double llik = 0;
   int nobs = w.length();
   int nlvar  = ncls.length();
   int nleaf = leaf.length();
   int nedge = edges.nrow();
   int nuniq = max(constr);

   List lrho(nuniq);
   List ltau(nedge);
   NumericVector lpi = pi_gnr(ncls[nlvar - 1], nobs) + log(w);

   List la(nlvar);
   List lb(nlvar);
   List ljbta(nedge);

   List lmpt(nlvar);
   List ljpst(nedge);


   std::vector<double*> prho(nleaf);
   std::vector<double*> ptau(nedge);
   double* ppi;

   std::vector<double*> pa(nlvar);
   std::vector<double*> pb(nlvar);
   std::vector<double*> pjbta(nedge);

   std::vector<double*> pmpt(nlvar);
   std::vector<double*> pjpst(nedge);

   for (int v = 0; v < nleaf; v ++) {
      NumericVector rho = rho_gnr(ncls[v], ncat[v]);
      prho[v] = rho.begin();
      lrho[v] = rho;
      NumericVector mpt(nobs * ncls[v]);
      pmpt[v] = mpt.begin();
      lmpt[v] = mpt;
   }

   for (int d = 0; d < nedge; d ++) {
      NumericVector tau = tau_gnr(ncls[edges(d, 0)], ncls[edges(d, 1)], nobs);
      ptau[d] = tau.begin();
      ltau[d] = tau;
      NumericVector jpst(nobs * ncls[edges(d, 0)] * ncls[edges(d, 1)]);
      pjpst[d] = jpst.begin();
      ljpst[d] = jpst;
      NumericVector jbta(nobs * ncls[edges(d, 1)]);
      pjbta[d] = jbta.begin();
      ljbta[d] = jbta;
   }

   for (int v = 0; v < nlvar; v ++) {
      NumericVector alpha(nobs * ncls[v]);
      NumericVector beta(nobs * ncls[v]);
      pa[v] = alpha.begin();
      pb[v] = beta.begin();
      la[v] = alpha;
      lb[v] = beta;
   }

   int iter = 0;
   double diff = R_PosInf;
   while (iter < max_iter || diff > tol) {
      iter ++;

      // (expectation-step)
      // initial beta, alpha
      int *py = y.begin();
      int *pnc = ncat.begin();
      for (int v = 0; v < nleaf; v ++) {
         msrPrd(py, prho[constr[v]], pb[v],
                ncls[v], nobs, nvar[v], pnc);
         py += nvar[v] * nobs;
         pnc += nvar[v];
      }

      // upward recursion
      for (int d = 0; d < nedge; d ++) {
         int u = edges(d, 0);
         int v = edges(d, 1);
         upRec(pb[v], pjbta[d], pb[u], ptau[d],
               nobs, ncls[u], ncls[v]);
      }

      // downward recursion
      la[nlvar - 1] = lpi;
      NumericVector beta1 = lb[nlvar - 1];
      NumericVector alpha1 = la[nlvar - 1];
      lmpt[nlvar - 1] = alpha1 + beta1;
      for (int d = nedge - 1; d == 0; d --) {
         int u = edges(d, 0);
         int v = edges(d, 1);
         int nl = ncls[v];  int nk = ncls[u];
         dnRec(pa[u], pa[v], pb[u], pb[v], pjbta[d],
               nobs, ncls[v], ncls[u], ptau[d],
               pmpt[u], pmpt[v], pjpst[d]);
      }

      List nlrho(nleaf);
      // (maximization-step)
      for (int v = 0; v < nleaf; v ++) {
         NumericVector opr = lrho[v];
         NumericVector npr(opr.length());
         int *py  = y.begin();
         double *rho = prho[constr[v]];
         double *pnpr = npr.begin();
         double *mpt = pmpt[v];
         for (int m = 0; m < nvar[v]; m ++) {
            for (int k = 0; k < ncls[v]; k ++) {
               for (int i = 0; i < nobs; i ++) {
                  if (py[i] > 0) prho[py[i]] += exp(mpt[v]);
                  else {
                     for (int r = 0; r < ncat[m]; r ++) {
                        pnpr[r] += exp(mpt[v] + prho[r]);
                     }
                  }
               }
               py += nobs;
            }
            pnpr += ncat[m];
            prho += ncat[m];
         }
      }
   }

   List post(nlvar);
   for (int v = 0; v < nlvar; v ++) {
      NumericVector alpha = la[v];
      NumericVector beta = lb[v];
      post[v] = alpha + beta;
   }

   List ret;
   ret["rho"] = lrho;
   ret["tau"] = ltau;
   ret["pi"] = lpi;
   ret["post"] = post;
   ret["llik"] = llik;

   return ret;
}


// beta is already included prevalence
void upRecS(
   double *nbeta, double *jbeta, double *beta,
   double *prv, double *lp,
   int nlvar, int nl, int *nk
) {
   double nu = 0;
   for (int j = 0; j < nlvar; j ++) {
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

   for (int l = 0; l < nlvar; l ++)
      nu += exp(nbeta[l]);
   for (int l = 0; l < nlvar; l ++)
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
