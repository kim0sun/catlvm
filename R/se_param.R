ginv <- function(x, tol = sqrt(.Machine$double.eps)) {
   xsvd <- svd(x)
   pos <- xsvd$d > tol * xsvd$d[1]

   if (all(pos))
      inv <- xsvd$v %*% (t(xsvd$u) / xsvd$d)
   else
      inv <- xsvd$v[, pos, drop = FALSE] %*%
         (t(xsvd$u[, pos, drop = FALSE]) / xsvd$d[pos])
   structure(inv, dimnames = dimnames(x)[2:1])
}

get_se <- function(vcov) {
   var <- diag(vcov)
   var[var < 0] = 0
   sqrt(var)
}

se_logit_fi <- function(scores, args) {
   score_list <- lapply(scores, do.call, what = rbind)
   score <- do.call(rbind, score_list)
   rownames(score) <- names(unlist(lapply(sapply(score_list, nrow), seq)))
   fi <- score %*% t(score)
   ginv(fi)
}

bdiag <- function(x) {
   nr <- sum(sapply(x, nrow))
   nc <- sum(sapply(x, ncol))
   ans <- matrix(0, nr, nc)
   ibegin <- 1; iend <- 0
   jbegin <- 1; jend <- 0

   for (m in x){
      iend <- iend + nrow(m)
      jend <- jend + ncol(m)
      ans[ibegin:iend, jbegin:jend] <- m
      ibegin <- ibegin + nrow(m)
      jbegin <- jbegin + ncol(m)
   }

   ans
}

split_tau <- function(x, args) {
   repeats <- rep(args$nclass_u, args$nclass_v)
   index <- rep(seq_len(sum(args$nclass_v)), unlist(repeats))
   split(unlist(x), index)
}

split_rho <- function(x, args) {
   repeats <- rep(args$ncat, args$nclass_leaf)
   index <- rep(seq(args$nvar %*% args$nclass_leaf), unlist(repeats))
   split(unlist(x), index)
}

jac_logistic <- function(x, simplify = TRUE) {
   ex <- exp(x)
   len <- length(ex)
   a <- diag(ex, len, len - 1)
   b <- (ex %o% ex)[, -len]
   a - b
}

se_transform <- function(log_par, logit_cov, args) {
   jacobians <- list()
   jacobians$pi <- bdiag(lapply(log_par$pi, jac_logistic))
   if (args$nlink > 0)
      jacobians$tau <- bdiag(lapply(split_tau(log_par$tau, args), jac_logistic))
   jacobians$rho <- bdiag(lapply(split_rho(log_par$rho, args), jac_logistic))

   jac <- bdiag(jacobians)
   rownames(jac) <- names(unlist(lapply(sapply(jacobians, nrow), seq)))
   covmat <- jac %*% logit_cov %*% t(jac)

   covmat
}
