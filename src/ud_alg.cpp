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
NumericVector ptr_gnr(int l, int k, int n) {
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

   return pv;
}

// [[Rcpp::export]]
NumericVector ptr_gnr2(int l, int k, int n) {
   NumericVector p(k);
   NumericVector pv(k * l);

   for (int i = 0; i < l; i ++) {
      p = runif(k, 0, 1);
      for (int j = 0; j < k; j ++) {
         pv[j + i * k] = log(p[j] / sum(p));
      }
   }

   return pv;
}


// [[Rcpp::export]]
NumericVector pem_gnr(int k, IntegerVector ncat) {
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
   int *y, double *lpre,
   int nk, int nobs,
   int nvar, int *ncat,
   double *beta
) {
   const double *lpre_init = lpre;
   for (int i = 0; i < nobs; i ++) {
      lpre = (double *)lpre_init;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0)
               beta[k] += lpre[y[m] - 1];

            lpre += ncat[m];
         }
      }
      y += nvar;
      beta += nk;
   }
}


// beta = [#class] * nobs
// jbeta = [#class * nlv] * nobs
// lbeta = [#cclass per nlv] * nobs
// lp = [#class * #pclass]
void upRec(
   double *beta, double *jbeta, std::vector<double*> lbeta,
   double *lprt, IntegerVector llv, int nlv, int nl, int *nk
) {
   for (int j = 0; j < nlv; j ++) {
      double *lbetaj = lbeta[llv[j]];
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk[j]; k ++)
            ml += exp(lprt[k] + lbetaj[k]);

         lprt += nk[j];
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
   }
}

void upRec2(
      double *beta, double *jbeta, double *lbeta,
      double *lprt, int nobs, int nl, int nk
) {
   for (int j = 0; j < nobs; j ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk; k ++)
            ml += exp(lprt[k] + lbeta[k]);

         lprt += nk;
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
      lbeta += nk;
   }
}

// nalpha = [#class] * nobs
// alpha = [#pclass] * nobs
// beta = [#pclass] * nobs
// ujbeta = [#pclass] * nobs
// lp = [#class * #pclass]
void dnRec(
   double *alpha, double *ualpha, double *ubeta, double *ujbeta, double *lprt,
   int nl, int nk
) {
   for (int k = 0; k < nk; k ++) {
      double val = 0;
      for (int l = 0; l < nl; l ++) {
         val += exp(lprt[k + l * nk] + ubeta[l] + ualpha[l] - ujbeta[l]);
      }
      alpha[k] = log(val);
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

List treeFit(
   IntegerVector y, IntegerVector w, int nobs,
   IntegerVector nvar, IntegerVector ncat, int nlv,
   IntegerMatrix edges, int nedge, IntegerVector ncls
) {
   List lpe(y.length());
   List lpt(nedge);

   List lalpa(nlv);
   List lbeta(nlv);
   List ljbta(nedge);

   List lmpt(y.length() + 1);
   List ljpt(nedge);


   std::vector<double*> ppe(y.length());
   std::vector<double*> ppt(nedge);
   double* ppp;

   std::vector<double*> palpa(nlv);
   std::vector<double*> pbeta(nlv);
   std::vector<double*> pjbta(nedge);

   std::vector<double*> pmpt(y.length() + 1);
   std::vector<double*> pjpt(nedge);

   for (int v = 0; v < y.length(); v ++) {
      NumericVector pem = pem_gnr(ncls[v], ncat[v]);
      ppe[v] = pem.begin();
      lpe[v] = pem;
      NumericVector mpt(nobs * ncls[v]);
      pmpt[v] = mpt.begin();
      lmpt[v] = mpt;
   }

   NumericVector mpt(nobs * ncls[nlv - 1]);
   pmpt[y.length()] = mpt.begin();
   lmpt[y.length()] = mpt;

   for (int d = 0; d < nedge; d ++) {
      NumericVector prt = ptr_gnr(ncls[edges(d, 0)], ncls[edges(d, 1)], nobs);
      ppt[d] = prt.begin();
      lpt[d] = prt;
      NumericVector jpt(nobs * ncls[edges(d, 0)] * ncls[edges(d, 1)]);
      pjpt[d] = prt.begin();
      ljpt[d] = prt;
      NumericVector jbta(nobs * ncls[edges(d, 1)]);
      pjbta[d] = jbta.begin();
      ljbta[d] = jbta;
   }

   for (int v = 0; v < nlv; v ++) {
      NumericVector alpha(nobs * ncls[v]);
      NumericVector beta(nobs * ncls[v]);
      palpa[v] = alpha.begin();
      pbeta[v] = beta.begin();
      lalpa[v] = alpha;
      lbeta[v] = beta;
   }


   int *py = y.begin();
   int *pnc = ncat.begin();
   for (int v = 0; v < y.length(); v ++) {
      msrPrd(py, ppe[v], ncls[v], nobs, nvar[v], pnc, pbeta[v]);
      py += nvar[v] * nobs;
      pnc += nvar[v];
   }

   for (int d = 0; d < nedge; d ++) {
      int vfrom = edges(d, 0); int vto = edges(d, 1);
      int nk = ncls[vfrom]; int nl = ncls[vto];
      double *lbeta = pbeta[vfrom]; double *beta = pbeta[vto];
      double *jbeta = pjbta[d];
      double *prt = ppt[d];

      for (int i = 0; i < nobs; i ++) {
         for (int l = 0; l < nl; l ++) {
            double ml = 0;
            for (int k = 0; k < nk; k ++)
               ml += exp(prt[k] + lbeta[k]);

            prt += nk;
            jbeta[l] = log(ml);
            beta[l] += log(ml);
         }
         jbeta += nl;
         beta  += nl;
      }
   }

   for (int d = nedge - 1; d == 0; d --) {
      int vfrom = edges(d, 1); int vto = edges(d, 0);
      int nl = ncls[vfrom];  int nk = ncls[vto];
      double *alpha = palpa[vto];
      double *ualpha = palpa[vfrom];
      double *beta = pbeta[vto];
      double *ubeta = pbeta[vfrom];
      double *jbeta = pjbta[d];
      double *mpt; double *jpt = pjpt[d];
      if (vto < y.length()) mpt = pmpt[vto];
      double *prt = ppt[d];

      for (int i = 0; i < nobs; i ++) {
         for (int k = 0; k < nk; k ++) {
            double jpr = 0;
            double mpr = 0;
            for (int l = 0; l < nl; l ++) {
               jpr = prt[k + l * nk] + ubeta[l] + ualpha[l] - jbeta[l];
               mpr += exp(jpr);
               jpt[k + l * nk] = jpr + beta[k];
            }
            alpha[k] = log(mpr);
            if (vto < y.length()) mpt[k] = alpha[k] + beta[k];
         }

         prt += nk * nl;
         jpt += nk * nl;
         if (vto < y.length()) mpt += nk;
         alpha += nk;
         ualpha += nl;
         beta += nk;
         ubeta += nl;
         jbeta += nl;
      }
   }

   List ret(1);
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
            vi[m] = j;
         }
         vi += k[i];
      }
   }
   return l;
}

// [[Rcpp::export]]
NumericVector test2(
      int n, IntegerVector k
) {
   NumericVector l(sum(k * k));
   IntegerVector index(n);

   int ind = 0;
   for (int i = 0; i < n; i ++) {
      index[i] = ind;
      ind += k[i] * k[i];
   }

   for (int i = 0; i < n; i ++) {
      for (int j = 0; j < k[i]; j ++) {
         for (int m = 0; m < k[i]; m ++) {
            l[index[i] + m + j * k[i]] = j;
         }
      }
   }

   return l;
}

// [[Rcpp::export]]
NumericVector t1(
      int I, int J, IntegerVector Jv, int niter
) {
   NumericVector v(J);
   double *vp = v.begin();

   for (int iter = 0; iter < niter; iter ++) {
      for (int i = 0; i < I; i ++) {
         for (int j = 0; j < J; j ++) {
            vp[j] += R::rnorm(0, 1);
         }
      }
   }

   return v;
}
// [[Rcpp::export]]
NumericVector t2(
      int I, int J, IntegerVector Jv, int niter
) {
   NumericVector v(J);
   double *vp = v.begin();

   for (int iter = 0; iter < niter; iter ++) {
      for (int i = 0; i < I; i ++) {
         for (int j = 0; j < Jv[i]; j ++) {
            vp[j] += R::rnorm(0, 1);
         }
      }
   }

   return v;
}
