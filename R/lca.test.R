library(glca)
data("gss08")

# MLCA
lca = glca(item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ 1,
           data = gss08, nclass = 3, n.init = 1, maxiter = 5)
par1 = lca$param$gamma
par2 = unlist(sapply(lca$param$rho[[1]], t))

pst = lcaf(t(as.matrix(lca$datalist$y[[1]])), rep(2, 6), 3, par1, par2)


# MLCA
lca = glca(item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ 1, DEGREE,
           data = gss08, nclass = 3, nclust = 2, n.init = 1, maxiter = 10)
par1 = lca$param$delta
par2 = c(t(lca$param$gamma))
par3 = sapply(lca$param$rho, t)

