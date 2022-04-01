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
NumericVector prb_gnr(int nk) {
   NumericVector p = runif(nk, 0, 1);
   return p / sum(p);
}

// [[Rcpp::export]]
NumericMatrix prb_gnrN(int nk, int N) {
   NumericVector p = runif(nk, 0, 1);
   NumericMatrix np(nk, N);

   for (int i = 0; i < N; i ++) {
      for (int j = 0; j < nk; j ++) {
         np(j, i) = p[j];
      }
   }
   return np;
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
   if (nvar == 0) return;
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

// List treeFit(
//    IntegerVector y, int nobs,
//    IntegerVector nvar, IntegerVector ncat,
//    int nlv, IntegerMatrix edges, int nedge,
//    IntegerVector ncls
// ) {
//    // List alpha(nlv);
//    // List beta(nlv);
//    // List jbeta(nedge);
//
//    List mpst(nlv);
//    List jpst(nedge);
//
//    std::vector<double*> plpre(y.length());
//    std::vector<double*> plprt(nedge);
//    double* plprp;
//
//    std::vector<double*> palpha(nlv);
//    std::vector<double*> pbeta(nlv);
//    std::vector<double*> pjbeta(nedge);
//
//    std::vector<double*> pmpst(y.length() + 1);
//    std::vector<double*> pjpst(nedge);
//
//
//    int *py = y.begin();
//    int *pnc = ncat.begin();
//    for (int v = 0; v < y.length(); v ++) {
//       msrPrd(py, plpre[v], ncls[v], nobs, nvar[v], pnc, pbeta[v]);
//       py += nvar[v] * nobs;
//       pnc += nvar[v];
//    }
//
//    for (int e = 0; e < nedge; e ++) {
//       int vfrom = edges(e, 0); int vto = edges(e, 1);
//       int nk = ncls[vfrom]; int nl = ncls[vto];
//       double *lbeta = pbeta[vfrom]; double *beta = pbeta[vto];
//       double *jbeta = pjbeta[e];
//       double *nlprt = plprt[e];
//
//       for (int i = 0; i < nobs; i ++) {
//          for (int l = 0; l < nl; l ++) {
//             double ml = 0;
//             for (int k = 0; k < nk; k ++)
//                ml += exp(nlprt[k] + lbeta[k]);
//
//             nlprt += nk;
//             jbeta[l] = log(ml);
//             beta[l] += log(ml);
//          }
//          jbeta += nl;
//          beta  += nl;
//       }
//    }
//
//    for (int e = nedge - 1; e == 0; e --) {
//       int vfrom = edges(e, 1); int vto = edges(e, 0);
//       int nl = ncls[vfrom];  int nk = ncls[vto];
//       double *alpha = palpha[vto];
//       double *ualpha = palpha[vfrom];
//       double *beta = pbeta[vto];
//       double *ubeta = pbeta[vfrom];
//       double *ujbeta = pjbeta[e];
//       double *mpst; double *jpst = pjpst[e];
//       if (vto < y.length()) mpst = pmpst[vto];
//       double *nlprt = plprt[e];
//
//       for (int i = 0; i < nobs; i ++) {
//          for (int k = 0; k < nk; k ++) {
//             double jpr = 0;
//             double mpr = 0;
//             for (int l = 0; l < nl; l ++) {
//                jpr = nlprt[k + l * nk] + ubeta[l] + ualpha[l] - ujbeta[l];
//                mpr += exp(jpr);
//                jpst[k + l * nk] = jpr + beta[k];
//             }
//             alpha[k] = log(mpr);
//             if (vto < y.length()) mpst[k] = alpha[k] + beta[k];
//          }
//
//          nlprt += nk * nl;
//          jpst += nk * nl;
//          if (vto < y.length()) mpst += nk;
//          alpha += nk;
//          ualpha += nl;
//          beta += nk;
//          ubeta += nl;
//          ujbeta += nl;
//       }
//    }
// }


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
NumericVector test(
   int I, int J, int niter
) {
   NumericVector v(I * J * niter);
   double *vp = v.begin();

   for (int iter = 0; iter < niter; iter ++) {
      for (int i = 0; i < I; i ++) {
         for (int j = 0; j < J; j ++) {
            vp[j] += R::rnorm(0, 1);
         }
         vp += J;
      }
   }

   return v;
}
// [[Rcpp::export]]
NumericVector test2(
      int I, int J, int niter
) {
   NumericVector v(I * J * niter);
   double *vp = v.begin();
   for (int iter = 0; iter < niter; iter ++) {
      for (int j = 0; j < J; j ++) {
         for (int i = 0; i < I; i ++) {
            vp[i * J + j] += R::rnorm(0, 1);
         }
      }
      vp += I * J;
   }
   return v;
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
