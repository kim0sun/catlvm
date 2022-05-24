output_param <- function(param, struct, args) {
   pi <- lapply(param$pi, exp)
   tau <- lapply(param$tau, exp)
   rho <- lapply(param$rho, exp)

   names(pi) <- struct$root

   for (d in seq_len(args$nedge)) {
      dimnames(tau[[d]]) <- list(
         seq(struct$nclass[args$u[d]]),
         seq(struct$nclass[args$v[d]])
      )
      names(dimnames(tau[[d]])) <- rev(unlist(struct$edges[d,]))
   }

   names(rho) <- letters[seq(args$nleaf_unique)]
   var <- split(struct$leaf, args$cstr_leaf)
   item <- split(struct$vars$manifest, args$cstr_leaf)
   for (v in seq_len(args$nleaf_unique)) {
      rho[[v]] <- matrix(rho[[v]], ncol = args$nclass_leaf[v])
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], seq_len),
         class = 1:ncol(rho[[v]])
      )

      var_index <- cumsum(c(1, args$ncat[[v]][-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1 (item ", seq(var_index), ")")
   }

   return(list(pi = pi, tau = tau, rho = rho))
}

output_logit <- function(lparam, struct, args) {
   pi <- lparam$pi
   tau <- lparam$tau
   rho <- lparam$rho

   for (r in seq_len(args$nroot)) {
      nclass <- args$nclass[args$root[r]]
      names(pi[[r]]) <- paste0(seq_len(nclass - 1), "/", nclass)
   }
   names(pi) <- struct$root

   for (d in seq_len(args$nedge)) {
      nk <- args$nclass[args$u[d]]
      nl <- args$nclass[args$v[d]]
      dimnames(tau[[d]]) <- list(
         paste0(seq(nk - 1), "/", nk), seq(nl)
      )
      names(dimnames(tau[[d]])) <- rev(unlist(struct$edges[d,]))
   }

   var <- split(struct$leaf, args$cstr_leaf)
   item <- split(struct$vars$manifest, args$cstr_leaf)
   for (v in seq(args$nleaf_unique)) {
      nclass <- args$nclass_leaf[v]
      rho[[v]] <- matrix(rho[[v]], ncol = nclass)
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], function(x)
            paste0(seq_len(x - 1), "/", x)),
         class = 1:ncol(rho[[v]])
      )

      var_index <- cumsum(c(1, (args$ncat[[v]] - 1)[-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1/", args$ncat[[v]], " (item ", seq(var_index), ")")
   }
   names(rho) <- letters[seq(args$nleaf_unique)]

   return(list(pi = pi, tau = tau, rho = rho))
}

output_posterior <- function(post, struct, data) {
   names(post) = struct$label
   lapply(post, function(x) {
      res <- exp(t(x))
      dimnames(res) <- list(data$dimnames[[1]], class = seq(ncol(res)))
      res
   })
}

