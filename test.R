source("~/Documents/Github/catlvm/R/formula.R")
source("~/Documents/Github/catlvm/R/args.R")
source("~/Documents/Github/catlvm/R/catlvm.R")
source("~/Documents/Github/catlvm/R/methods.R")
source("~/Documents/Github/catlvm/R/simulate.catlvm.R")
source("~/Documents/Github/catlvm/R/data.R")
source("~/Documents/Github/catlvm/R/plot.catlvm.R")
source("~/Documents/Github/catlvm/R/estimate.R")
source("~/Documents/Github/catlvm/R/param_validate.R")
source("~/Documents/Github/catlvm/R/output_fit.R")
source("~/Documents/Github/catlvm/R/transform_param.R")
source("~/Documents/Github/catlvm/R/catlvm.control.R")
Rcpp::sourceCpp("~/Documents/Github/catlvm/src/ud_alg.cpp")
library(dplyr)

lca   = catlvm(L1[2] ~ X1 + X2 + X3)
lcas  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[2] ~ Y1 + Y2 + Y3)
jlca  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[4] ~ Y1 + Y2 + Y3,
               L3[4] ~ Z1 + Z2 + Z3,
               JC[3] ~ L1 + L2 + L3)
lcpa  = catlvm(L1[2] ~ X1 + X2 + X3,
               L2[2] ~ Y1 + Y2 + Y3,
               L3[2] ~ Z1 + Z2 + Z3,
               PF[2] ~ L1 + L2 + L3,
               constraints = list(c("L1", "L2", "L3"),
                                  c("PF~L1", "PF->L2", "PF ->L3")))
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
jjcpa = catlvm(L1[2] ~ X11 + X21 + X31,
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
               JP[2] ~ J1 + J2 + J3)
lta = catlvm(L1[3] ~ X11 + X21 + X31,
             L2[3] ~ X12 + X22 + X32,
             L3[3] ~ X13 + X23 + X33,
             L1 ~ L2, L2 ~ L3,
             constraints = list(c("L1", "L2", "L3")))
jlta = catlvm(L1[3] ~ X11 + X21 + X31,
              M1[3] ~ Y11 + Y21 + Y31,
              N1[3] ~ Z11 + Z21 + Z31,
              L2[3] ~ X12 + X22 + X32,
              M2[3] ~ Y12 + Y22 + Y32,
              N2[3] ~ Z12 + Z22 + Z32,
              L3[3] ~ X13 + X23 + X33,
              M3[3] ~ Y13 + Y23 + Y33,
              N3[3] ~ Z13 + Z23 + Z33,
              J1[3] ~ L1 + M1 + N1,
              J2[3] ~ L2 + M2 + N2,
              J3[3] ~ L3 + M3 + N3,
              J1 ~ J2, J2 ~ J3,
             constraints = list(c("L1", "L2", "L3"),
                                c("M1", "M2", "N3"),
                                c("N1", "N2", "N3")))
# constraint to 0
# constraint edge X
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

library(dplyr)
object <- jlcpa
sim <- object %>% simulate(1000)
fit <- object %>% estimate(data = sim$response, method = "hybrid")
fit$fit$posterior
fit$fit$estimates

library(randomLCA)
data(symptoms)
dat = symptoms[rep(seq(symptoms$Freq), symptoms$Freq),]


lcpa_sym <- catlvm(
   LC1[3] ~ Nightcough.13 + Wheeze.13 + Itchyrash.13 + FlexDerma.13,
   LC2[3] ~ Nightcough.45 + Wheeze.45 + Itchyrash.45 + FlexDerma.45,
   LC3[3] ~ Nightcough.6  + Wheeze.6  + Itchyrash.6  + FlexDerma.6,
   LC4[3] ~ Nightcough.7  + Wheeze.7  + Itchyrash.7  + FlexDerma.7,
   LCP[2] ~ LC1 + LC2 + LC3 + LC4,
   constraints = list(c("LC1", "LC2", "LC3", "LC4"))
)

lcpa_sym2 <- catlvm(
   LC1[3] ~ Nightcough.13 + Wheeze.13 + Itchyrash.13 + FlexDerma.13,
   LC2[3] ~ Nightcough.45 + Wheeze.45 + Itchyrash.45 + FlexDerma.45,
   LC3[3] ~ Nightcough.6  + Wheeze.6  + Itchyrash.6  + FlexDerma.6,
   LC4[3] ~ Nightcough.7  + Wheeze.7  + Itchyrash.7  + FlexDerma.7,
   LCP[2] ~ LC1 + LC2 + LC3 + LC4,
   constraints = list(c("LC1", "LC2", "LC3", "LC4")),
   data = dat
)

object = lcpa_sym %>% estimate(data = dat)
predict(object, "LC1")
a = rho(object)

lcpa_sym %>% estimate(data = dat)
undebug(estimate.catlvm)

library(randomLCA)
data(symptoms)
