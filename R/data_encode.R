rep_row <- function(x, lev) {
   sb <- as.matrix(expand.grid(lev[x == 0]))
   m <- t(replicate(nrow(sb), x))
   m[, x == 0] <- sb
   m
}

split_by_name <- function(x, f) lapply(f, function(i) x[i])

stretch_data <- function(items, mf) {
   m <- split_by_name(mf, items)
   y <- lapply(m, function(x) t(sapply(x, as.numeric)))
   unlist(unname(y))
}

proc_saturated <- function(mf, ncat) {
   lev <- lapply(ncat, seq_len)
   if (all(unlist(mf) > 0)) {
      yobs <- mf
      y_aggr <- aggregate(numeric(nrow(yobs)), yobs, length)
      y_unique <- y_aggr[, -ncol(y_aggr)]
      freq <- y_aggr[, ncol(y_aggr)]
      loglik <- sum(freq * log(freq / sum(freq)))

      res <- list(y = y_unique, freq = freq, loglik = loglik)
   } else {
      na_ind <- rowSums(mf == 0) > 0
      yobs <- mf[!na_ind, , drop = FALSE]
      ymis <- mf[ na_ind, , drop = FALSE]

      yobs0 <- aggregate(numeric(nrow(yobs)), yobs, length)
      ymis0 <- aggregate(numeric(nrow(ymis)), ymis, length)
      y_expand <- do.call(rbind, apply(
         ymis0[, -ncol(ymis0)], 1, rep_row, lev, simplify = FALSE))
      y0 <- rbind(yobs0, cbind(y_expand, x = 0))

      y_aggr <- aggregate(y0[[ncol(y0)]], y0[-ncol(y0)], sum)
      y_unique <- y_aggr[, -ncol(y0)]
      freq <- y_aggr[, ncol(y0)]

      nrep <- apply(ymis0[, -ncol(ymis0)], 1,
                    function(x) prod(ncat[x == 0]))
      miss <- match(apply(y_expand, 1, paste, collapse = ""),
                    apply(y_unique, 1, paste, collapse = "")) - 1

      misEM <- calcfreq(miss, nrep, nrow(ymis0),
                        ymis0[, ncol(ymis0)], freq,
                        nrow(y_aggr), nrow(mf), 1e-5, 100)

      loglik <- sum(misEM$freq * log(misEM$theta)) + misEM$loglik

      res <- list(y = y_unique, freq = calc_mis$freq, loglik = loglik)
   }

   res
}

proc_data <- function(data, model) {
   if (is.null(data)) return(NULL)
   items <- model$vars$manifest
   f <- paste("~", paste(unlist(items), collapse = "+"))
   mf <- model.frame(formula(f), data)
   mf[] <- lapply(mf, factor)
   level <- lapply(mf, levels)
   ncat <- sapply(mf, nlevels)
   mf[] <- lapply(mf, as.numeric)
   mf[is.na(mf)] = 0
   dims <- dimnames(mf)

   if (any(!unlist(items) %in% dims[[2]]))
      stop("Following manifest variables not found:\n ",
           paste(unlist(items)[!unlist(items) %in% dims[[2]]],
                 collapse = " "))

   saturated <- proc_saturated(mf, ncat)
   saturated$y <- stretch_data(items, saturated$y)

   list(y = stretch_data(items, mf),
        nobs = nrow(mf),
        ncat = split_by_name(ncat, items),
        level = split_by_name(level, items),
        saturated = saturated,
        dimnames = dims)
}

