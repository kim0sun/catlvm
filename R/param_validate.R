pi_valid <- function(pi, nclass, gnr = FALSE) {
   if (is.numeric(pi)) {
      lenTF <- length(pi) == c(nclass)
      sumTF <- abs(sum(pi) - 1) < .Machine$double.eps^0.5
      if (lenTF && sumTF)
         return(if (gnr) log(pi) else TRUE)
      else return(if (gnr) pi_gnr(nclass)
                  else FALSE)
   } else return(if (gnr) pi_gnr(nclass)
                 else FALSE)
}

tau_valid <- function(tau, nk, nl, gnr = FALSE) {
   if (is.matrix(tau)) {
      dimTF <- all(dim(tau) == c(nk, nl))
      sumTF <- all(abs(colSums(tau) - 1) < .Machine$double.eps^0.5)
      if (dimTF && sumTF)
         return(if (gnr) log(tau) else TRUE)
      else return(if (gnr) tau_gnr(nk, nl) else FALSE)
   } else return(if (gnr) tau_gnr(nk, nl) else FALSE)
}

rho_valid <- function(rho, nclass, ncat, gnr = FALSE) {
   if (is.list(rho)) {
      dim1 <- lapply(rho, dim)
      dim2 <- lapply(ncat, function(m) c(m, nclass))
      dimTF <- identical(dim1, dim2)
      sumTF <- all(sapply(rho, function(x)
         all(abs(colSums(x) - 1) < .Machine$double.eps^0.5)))
      if (dimTF && sumTF)
         return(if (gnr) log(unlist(rho)) else TRUE)
      else return(if (gnr) rho_gnr(nclass, ncat) else FALSE)
   } else if (is.numeric(rho)) {
      lrho <- matrix(rho, ncol = nclass)
      ind <- rep(seq(length(ncat)), ncat)
      lenTF <- length(rho) == nclass * sum(ncat)
      sumTF <- all(apply(lrho, 2, function(x)
         all(abs(sapply(split(x, ind), sum) - 1) < .Machine$double.eps^0.5)))

      if (lenTF && sumTF) return(if (gnr) log(rho) else TRUE)
      else return(if (gnr) rho_gnr(nclass, ncat) else FALSE)
   } else return(if (gnr) rho_gnr(nclass, ncat) else FALSE)
}

init_validate <- function(init.param, args) {
   valid <- c()
   pi <- init.param$pi
   tau <- init.param$tau
   rho <- init.param$rho

   if (length(pi) < args$nroot) valid["pi"] <- FALSE
   else valid["pi"] <- all(sapply(seq_len(args$nroot), function(r)
      pi_valid(pi[[r]], args$nclass[args$root[r]])))

   if (length(tau) < args$nedge) valid["tau"] <- FALSE
   else valid["tau"] <- all(sapply(seq_len(args$nedge), function(d)
      tau_valid(tau[[d]], args$nclass[args$u[d]], args$nclass[args$v[d]])))

   if (length(rho) < args$nleaf_unique) valid["rho"] <- FALSE
   else valid["rho"] <- all(sapply(seq_len(args$nleaf_unique), function(v)
      rho_valid(rho[[v]], args$nclass_leaf[v], args$ncat[[v]])))

   valid
}
