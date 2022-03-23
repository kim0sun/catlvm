#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


// cll (ncls * nobs/nlv)
//

double lnadd(double lx, double ly) {
   if (lx == R_NegInf) return ly;
   if (ly == R_NegInf) return lx;

   double lxy = 0;
   if (lx < ly)
      lxy = ly + log(1 + exp(lx - ly));
   else
      lxy = lx + log(1 + exp(ly - lx));

   return lxy;
}

void msr_obs(
   int *y, double *lp,
   int ncls, int nobs,
   int nvar, int *ncat,
   double *cll
) {
   const double *lp_init = lp;
   for (int i = 0; i < nobs; i ++) {
      lp = (double *)lp_init;
      for (int m = 0; m < nvar; m ++) {
         for (int k = 0; k < ncls; k ++) {
            if (y[m] > 0)
               cll[k] += lp[y[m] - 1];
            if (i == 0)
               Rcout << ncat[m] << std::endl;
            lp += ncat[m];
         }
      }
      y   += nvar;
      cll += ncls;
   }
}

void msr_lv(
   double *clls, double *lp,
   int nlv, int *nobs,
   int ncls, int *nclss,
   double *dpr, double *cll
) {
   for (int j = 0; j < nlv; j ++) {
      const double *lp_init = lp;
      for (int i = 0; i < nobs[j]; i ++) {
         lp = (double *)lp_init;
         for (int l = 0; l < ncls; l ++) {
            double mll = R_NegInf;
            for (int k = 0; k < nclss[j]; k ++)
               mll = lnadd(mll, lp[k] + clls[k]);
            for (int k = 0; k < nclss[j]; k ++)
               dpr[k] = lp[k] + clls[k] - mll;

            lp   += nclss[j];
            dpr  += nclss[j];
            cll[l] += mll;
         }
         clls += nclss[j];
      }
      cll += ncls;
   }
}

double prd_prev(
   double *cll, double *lp,
   int ncls, int nobs, double *pst
) {
   double ll = 0;
   for (int i = 0; i < nobs; i ++) {
      double mll = R_NegInf;
      for (int k = 0; k < ncls; k ++) {
         mll = lnadd(mll, lp[k] + cll[k]);
      }
      for (int k = 0; k < ncls; k ++) {
         pst[k] = lp[k] + cll[k] - mll;
      }
      cll += ncls;
      pst += ncls;
      ll  += mll;
   }
   return ll;
}

// [[Rcpp::export]]
NumericMatrix lcaf(
   IntegerMatrix y, IntegerVector ncat, int ncls,
   NumericVector prev, NumericVector par
) {
   NumericMatrix cll(ncls, y.ncol());
   NumericMatrix pst(ncls, y.ncol());
   msr_obs(y.begin(), par.begin(), ncls,
           y.ncol(), y.nrow(), ncat.begin(),
           cll.begin());
   double ll = prd_prev(cll.begin(), prev.begin(),
            ncls, y.ncol(), pst.begin());
   Rcout << ll << std::endl;
   return pst;
}


// [[Rcpp::export]]
NumericMatrix mlcaf(
      IntegerMatrix y, IntegerVector ncat, int ncls,
      NumericVector prev, NumericVector par
) {
   NumericMatrix cll(ncls, y.ncol());
   NumericMatrix pst(ncls, y.ncol());
   msr_obs(y.begin(), par.begin(), ncls,
           y.ncol(), y.nrow(), ncat.begin(),
           cll.begin());
   double ll = prd_prev(cll.begin(), prev.begin(),
                        ncls, y.ncol(), pst.begin());
   Rcout << ll << std::endl;
   return pst;
}

bool a(){
   return(true);
}

