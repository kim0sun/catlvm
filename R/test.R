source("~/Documents/Github/catlvm/R/formula.R")
source("~/Documents/Github/catlvm/R/constraints.R")
source("~/Documents/Github/catlvm/R/args.R")
source("~/Documents/Github/catlvm/R/catlvm.R")
source("~/Documents/Github/catlvm/R/Methods.R")
Rcpp::sourceCpp("~/Documents/Github/catlvm/src/ud_alg2.cpp")

lca   = catlvm(L1[2] ~ X1 + X2 + X3)
lcas  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[2] ~ Y1 + Y2 + Y3,
               L3[2] ~ Z1 + Z2 + Z3,
               L4[2] ~ W1 + W2 + W3)
jlca  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[2] ~ Y1 + Y2 + Y3,
               L3[2] ~ Z1 + Z2 + Z3,
               JC[2] ~ L1 + L2 + L3)
lcpa  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[2] ~ Y1 + Y2 + Y3,
               L3[2] ~ Z1 + Z2 + Z3,
               PF[2] ~ L1 + L2 + L3,
               constraints = list(c("L1", "L2", "L3")))
jlcpa = catlvm(L1[2] ~ X11 + X21 + X31,
               M1[2] ~ Y11 + Y21 + Y31,
               N1[2] ~ Z11 + Z21 + Z31,
               L2[2] ~ X12 + X22 + X32,
               M2[2] ~ Y12 + Y22 + Y32,
               N2[2] ~ Z12 + Z22 + Z32,
               L3[2] ~ X13 + X23 + X33,
               M3[2] ~ Y13 + Y23 + Y33,
               N3[2] ~ Z13 + Z23 + Z33,
               J1[2] ~ L1 + M1 + N1,
               J2[2] ~ L2 + M2 + N2,
               J3[2] ~ L3 + M3 + N3,
               JP[2] ~ J1 + J2 + J3,
               constraints = list(c("L1", "L2", "L3"),
                                  c("M1", "M2", "M3"),
                                  c("N1", "N2", "N3")))
lta = catlvm(L1[2] ~ X11 + X21 + X31,
             L2[2] ~ X12 + X22 + X32,
             L3[2] ~ X13 + X23 + X33,
             L1 ~ L2, L2 ~ L3,
             constraints = list(c("L1", "L2", "L3")))
lcawg = catlvm(LG[2] ~ Z1 + Z2 + Z3,
               LC[2] ~ X1 + X2 + X3,
               LG ~ LC)
lcpawg = catlvm(LG[2] ~ Z1 + Z2 + Z3,
                LG ~ P1,
                L1[2] ~ X11 + X12 + X13,
                L2[2] ~ X21 + X22 + X23,
                L3[2] ~ X31 + X32 + X33,
                P1[2] ~ L1 + L2 + L3,
                constraints = list(c("L1", "L2", "L3")))

lca; plot(lca)
lcpa; plot(lcpa)
jlca; plot(jlca, abbreviation = TRUE)
jlcpa; plot(jlcpa, abbreviation = TRUE)
lta; plot(lta)

lta %>% estimate(data = response)
lta %>% estimate(data = response) %>%
   regression(L1 ~ Cov1 + Cov2,
              L2 ~ Cov1 + Cov2,
              L3 ~ Cov1 + Cov2,
              data = response)

nsim = 1000
niter = 200
{
   object = lta
   y = simulate(object, nsim)
   par = list()
   for (iter in 1:niter) {
      cat(iter, "\r")
      y = simulate(object, nsim,
                   tau = lapply(y$params$tau, exp),
                   rho = lapply(y$params$rho, exp))
      fit = emFit(y = unlist(y$response),
                  nobs = y$args$nobs, nvar = y$args$nvar, ncat = y$args$ncat,
                  nlv = y$args$nlv, root = y$args$root - 1, leaf =  y$args$leaf - 1,
                  cstr_leaf = y$args$cstr_leaf - 1,
                  ulv = y$args$u - 1, vlv = y$args$v - 1,
                  cstr_root = y$args$cstr_root - 1, cstr_edge = y$args$cstr_edge - 1,
                  nclass = y$args$nclass, nclass_leaf = y$args$nclass_leaf,
                  nclass_u = y$args$nclass_u, nclass_v = y$args$nclass_v,
                  init = rep(TRUE, 3), init_param =  y$params,
                  max_iter = 1e4, tol = 1e-3, verbose = FALSE)
      par[[iter]] = fit$params
   }
}

lapply(y$params$pi, exp)
lapply(1:length(y$params$pi), function(i) Reduce("+", lapply(par, function(x) exp(x$pi[[i]]))) / niter)
lapply(y$params$tau, exp)
lapply(1:length(y$params$tau), function(i) Reduce("+", lapply(par, function(x) exp(x$tau[[i]]))) / niter)
lapply(y$params$rho, exp)
lapply(1:length(y$params$rho), function(i) Reduce("+", lapply(par, function(x) exp(x$rho[[i]]))) / niter)
#
# lapply(y$params$pi, function(x) round(exp(x), 3))
# lapply(fit$param$pi, function(x) round(exp(x), 3))
# lapply(y$params$tau, function(x) round(exp(x), 3))
# lapply(fit$param$tau, function(x) round(exp(x), 3))
# lapply(y$params$rho, function(x) round(exp(x), 3))
# lapply(fit$param$rho, function(x) round(exp(x), 3))

