#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

class catlv{
public:
   int depth;
   bool isleaf;
   bool isroot;

   int parent;
   int *children;

   double *beta;
   double *alpha;

   catlv(List _info) {
      depth = _info["depth"];
      isleaf = _info["isleaf"];
      isroot = _info["isroot"];
      parent = _info["parent"];
      children = _info["children"];
   }

   void upRec();
   void dnRec();
   void udPrd();
   void upRecS();
   void dnRecS();
   void udPrdS();
};

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
   int *y, double *lp,
   int nk, int nobs,
   int nvar, int *ncat,
   double *beta
) {
   const double *lp_init = lp;
   for (int i = 0; i < nobs; i ++) {
      lp = (double *)lp_init;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < nk; k ++) {
            if (y[m] > 0)
               beta[k] += lp[y[m] - 1];

            lp += ncat[m];
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
   double *beta, double *jbeta, double *lbeta, double *lp,
   int nlv, int nl, int *nk
) {
   for (int j = 0; j < nlv; j ++) {
      for (int l = 0; l < nl; l ++) {
         double ml = 0;
         for (int k = 0; k < nk[j]; k ++)
            ml += exp(lp[k] + lbeta[k]);

         lp += nk[j];
         jbeta[l] = log(ml);
         beta[l] += log(ml);
      }
      jbeta += nl;
      lbeta += nk[j];
   }
}


// nalpha = [#class] * nobs
// alpha = [#pclass] * nobs
// beta = [#pclass] * nobs
// ujbeta = [#pclass] * nobs
// lp = [#class * #pclass]
void dnRec(
   double *alpha, double *ualpha, double *ubeta, double *ujbeta, double *lp,
   int nl, int nk
) {
   for (int k = 0; k < nk; k ++) {
      double val = 0;
      for (int l = 0; l < nl; l ++) {
         val += exp(lp[k + l * nk] + ubeta[l] + ualpha[l] - ujbeta[l]);
      }
      alpha[k] = log(val);
   }
}

void udPrd(
   double *mar, double *joint,
   double *alpha, double *ualpha,
   double *beta, double *jbeta, double *ubeta,
   double *lp, int nl, int nk
) {
   for (int l = 0; l < nl; l ++) {
      for (int k = 0; k < nk; k ++) {
         joint[l] = beta[k] + lp[k] + ualpha[l] + ubeta[l] - jbeta[l];

         if (l == 0) {
            mar[k] = alpha[k] + beta[k];
         }
      }
      lp += nk;
      joint += nk;
   }
}


void treePst(
   int ndepth,
   IntegerVector nk,
   IntegerVector plv,
   List clv,
   List beta
) {
   // measure leaf lv

   // elevate internal lv -> root lv

   // decline internal lv -> leaf lv

   // depth declining

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
