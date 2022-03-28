library(glca)
data("gss08")

# MLCA
nc = 5
lca = glca(item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ 1,
           data = gss08, nclass = nc, n.init = 1)
par1 = log(lca$param$gamma)
par2 = log(unlist(sapply(lca$param$rho[[1]], t)))

yt = t(as.matrix(lca$datalist$y[[1]]))

pst = lcaf(yt, rep(2, 6), nc, par1, par2)
round(t(exp(pst))[1:10,], 3)
round(lca$posterior[[1]][1:10,], 3)
lca$gof$loglik

# MLCA
mlca = glca(item(DEFECT, HLTH, RAPE, POOR, SINGLE, NOMORE) ~ 1, DEGREE,
           data = gss08, nclass = 3, nclust = 2, n.init = 1)
par1 = log(mlca$param$delta)
par2 = log(c(t(mlca$param$gamma)))
par3 = log(sapply(mlca$param$rho, t))

yt = lapply(mlca$datalist$y, function(y) t(as.matrix(y)))
gp = mlca$datalist$group

lcaf(yt[[1]], rep(2, 6), 3, par)
