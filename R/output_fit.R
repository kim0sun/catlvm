output_param <- function(struct, args, params) {
   pi <- lapply(params$pi, exp)
   tau <- lapply(params$tau, exp)
   rho <- lapply(params$rho, exp)

   names(pi) <- struct$root

   for (d in seq_len(length(tau))) {
      dimnames(tau[[d]]) <- list(
         seq(struct$nclass[args$u[d]]),
         seq(struct$nclass[args$v[d]])
      )
      names(dimnames(tau[[d]])) <- rev(unlist(struct$edges[d,]))
      class(tau[[d]]) = "catlvm.tau"
   }

   names(rho) <- letters[seq(args$nleaf_unique)]
   var <- split(struct$leaf, args$cstr_leaf)
   item <- split(struct$vars$manifest, args$cstr_leaf)
   for (v in seq(args$nleaf_unique)) {
      rho[[v]] <- matrix(rho[[v]], ncol = args$nclass_leaf[v])
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], seq_len),
         class = 1:ncol(rho[[v]])
      )

      attr(rho[[v]], "variables") <- as.character(var[[v]])
      attr(rho[[v]], "items") <- item[[v]]

      var_index <- cumsum(c(1, args$ncat[[v]][-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1 (item ", seq(var_index), ")")

      class(rho[[v]]) = "catlvm.irp"
   }

   return(list(pi = pi, tau = tau, rho = rho))
}

output_logit <- function(struct, args, params) {
   pi <- lapply(params$pi, exp)
   tau <- lapply(params$tau, exp)
   rho <- lapply(params$rho, exp)

   names(pi) <- struct$root

   for (d in seq_len(length(tau))) {
      dimnames(tau[[d]]) <- list(
         seq(struct$nclass[args$u[d]]),
         seq(struct$nclass[args$v[d]])
      )
      names(dimnames(tau[[d]])) <- rev(unlist(struct$edges[d,]))
      class(tau[[d]]) = "catlvm.tau"
   }

   names(rho) <- letters[seq(args$nleaf_unique)]
   var <- split(struct$leaf, args$cstr_leaf)
   item <- split(struct$vars$manifest, args$cstr_leaf)
   for (v in seq(args$nleaf_unique)) {
      rho[[v]] <- matrix(rho[[v]], ncol = args$nclass_leaf[v])
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], seq_len),
         class = 1:ncol(rho[[v]])
      )

      attr(rho[[v]], "variables") <- as.character(var[[v]])
      attr(rho[[v]], "items") <- item[[v]]

      var_index <- cumsum(c(1, args$ncat[[v]][-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1 (item ", seq(var_index), ")")

      class(rho[[v]]) = "catlvm.irp"
   }

   return(list(pi = pi, tau = tau, rho = rho))
}

output_se <- function(struct, args, params) {
   pi <- lapply(params$pi, exp)
   tau <- lapply(params$tau, exp)
   rho <- lapply(params$rho, exp)

   names(pi) <- struct$root

   for (d in seq_len(length(tau))) {
      dimnames(tau[[d]]) <- list(
         seq(struct$nclass[args$u[d]]),
         seq(struct$nclass[args$v[d]])
      )
      names(dimnames(tau[[d]])) <- rev(unlist(struct$edges[d,]))
      class(tau[[d]]) = "catlvm.tau"
   }

   names(rho) <- letters[seq(args$nleaf_unique)]
   var <- split(struct$leaf, args$cstr_leaf)
   item <- split(struct$vars$manifest, args$cstr_leaf)
   for (v in seq(args$nleaf_unique)) {
      rho[[v]] <- matrix(rho[[v]], ncol = args$nclass_leaf[v])
      dimnames(rho[[v]]) <- list(
         reponse = sapply(args$ncat[[v]], seq_len),
         class = 1:ncol(rho[[v]])
      )

      attr(rho[[v]], "variables") <- as.character(var[[v]])
      attr(rho[[v]], "items") <- item[[v]]

      var_index <- cumsum(c(1, args$ncat[[v]][-args$nvar[v]]))
      rownames(rho[[v]])[var_index] <-
         paste0("1 (item ", seq(var_index), ")")

      class(rho[[v]]) = "catlvm.irp"
   }

   return(list(pi = pi, tau = tau, rho = rho))
}

output_posterior <- function(struct, args, data, post) {
   names(post) = struct$label
   lapply(post, function(x) {
      res <- exp(t(x))
      dimnames(res) <- list(data$dimnames[[1]], class = seq(ncol(res)))
      res
   })
}

