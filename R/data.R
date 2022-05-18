rep_row <- function(x, lev) {
   sb <- as.matrix(expand.grid(lev[x == 0]))
   m <- t(replicate(nrow(sb), x))
   m[, x == 0] <- sb

   cbind(m, x = 0)
}

match_na <- function(x, y) {
   px <- paste(x[x > 0], collapse = "")
   py <- apply(y[, x > 0], 1, paste, collapse = "")
   which(px == py)
}

split_by_name <- function(x, f) lapply(f, function(i) x[i])

stretch_data <- function(items, mf) {
   m <- split_by_name(mf, items)
   y <- lapply(m, function(x) t(sapply(x, as.numeric)))
   unlist(unname(y))
}

proc_saturated <- function(mf, items, ncat) {
   lev <- lapply(ncat, seq_len)
   if (all(unlist(mf) > 0)) {
      yobs <- mf
      yn <- aggregate(numeric(nrow(yobs)), yobs, length)
      unique_y <- yn[, -ncol(yn)]
      freq <- yn[, ncol(yn)]
      loglik <- sum(freq * log(freq / sum(freq)))

      res <- list(y = stretch_data(items, unique_y),
                  nobs = nrow(unique_y),
                  freq = freq, loglik = loglik)
   } else {
      na_ind <- rowSums(mf == 0) > 0
      yobs <- mf[!na_ind, , drop = FALSE]
      ymis <- mf[ na_ind, , drop = FALSE]

      yobs0 <- aggregate(numeric(nrow(yobs)), yobs, length)
      ymis0 <- aggregate(numeric(nrow(ymis)), ymis, length)
      expand_y <- do.call(rbind, apply(ymis0[, -ncol(ymis0)], 1, rep_row,
                                       lev, simplify = FALSE))
      y0 <- rbind(yobs0, expand_y)

      yn <- aggregate(y0[[ncol(y0)]], y0[-ncol(y0)], sum)
      unique_y <- yn[, -ncol(y0)]
      freq <- yn[, ncol(y0)]

      mis_patt <- apply(ymis0[, -ncol(ymis0)], 1, match_na,
                        unique_y, simplify = FALSE)
      nrep <- as.numeric(sapply(mis_patt, length))
      miss <- unlist(mis_patt) - 1

      calc_mis <- calcfreq(miss, nrep, nrow(ymis0),
                           ymis0[, ncol(ymis0)], freq,
                           nrow(yn), nrow(mf), 1e-5, 100)

      theta <- calc_mis$freq / sum(calc_mis$freq)
      loglik <- sum(freq * log(theta)) + calc_mis$loglik

      res <- list(y = stretch_data(items, unique_y),
                  nobs = nrow(unique_y),
                  freq = calc_mis$freq,
                  loglik = loglik)
   }

   res
}

proc_data <- function(data, struct) {
   if (is.null(data)) return(NULL)
   items <- struct$vars$manifest
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

   saturated_model <- proc_saturated(mf, items, ncat)

   list(y = stretch_data(items, mf),
        nobs = nrow(mf),
        ncat = split_by_name(ncat, items),
        level = split_by_name(level, items),
        saturated = saturated_model,
        dimnames = dims)
}

