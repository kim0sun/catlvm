rep_row <- function(x, level) {
   na <- which(is.na(x))
   replace <- as.matrix(expand.grid(level[na]))
   m <- t(replicate(nrow(replace), x))
   m[, na] <- replace
   cbind(m, x = 0)
}

match_na <- function(x, y) {
   na <- is.na(x)
   px <- paste(x[!na], collapse = "")
   py <- apply(y[, !na], 1, paste, collapse = "")
   which(px == py)
}

split_by_name <- function(x, f) lapply(f, function(i) x[i])

stretch_data <- function(items, mf) {
   m <- split_by_name(mf, items)
   y <- lapply(m, function(x) t(sapply(x, as.numeric)))
   unlist(unname(y))
}

proc_saturated <- function(mf, items) {
   mf[] <- lapply(mf, as.numeric)
   lev <- lapply(mf, function(x) seq_len(max(x)))

   if (anyNA(mf)) {
      na_ind <- rowSums(is.na(mf)) > 0
      yobs <- mf[!na_ind, , drop = FALSE]
      ymis <- mf[ na_ind, , drop = FALSE]
   } else {
      yobs <- mf
      ymis <- mf[0,]
   }

   yobs0 <- aggregate(numeric(nrow(yobs)), yobs, length)
   ymis0 <- do.call(rbind, apply(ymis, 1, rep_row, lev, simplify = FALSE))
   y0 <- rbind(yobs0, ymis0)
   aggr_y <- aggregate(y0[[ncol(y0)]], y0[-ncol(y0)], sum)
   uniq_y <- aggr_y[, -ncol(y0)]

   mis_patt <- apply(ymis, 1, match_na, uniq_y)
   nrep <- as.numeric(sapply(mis_patt, length))
   miss <- unlist(mis_patt) - 1

   freq <- calcfreq(miss, nrep, nrow(ymis), aggr_y[[ncol(aggr_y)]],
                    nrow(aggr_y), nrow(mf), 1e-5, 100)
   term <- terms_data(mi, lapply(uniq_y, factor))

   list(y = stretch_data(items, uniq_y),
        nobs = nrow(uniq_y),
        freq = freq)
}

proc_data <- function(data, struct) {
   if (is.null(data)) return(NULL)
   items <- struct$vars$manifest
   f <- paste("~", paste(unlist(items), collapse = "+"))
   mf <- model.frame(formula(f), data)
   mf[] <- lapply(mf, factor)
   level <- lapply(mf, levels)
   ncat <- sapply(mf, nlevels)
   dims <- dimnames(mf)

   if (any(!unlist(items) %in% dims[[2]]))
      stop("Following manifest variables not found:\n ",
           paste(unlist(items)[!unlist(items) %in% dims[[2]]],
                 collapse = " "))

   saturated_model <- proc_saturated(mf, items)

   list(y = stretch_data(items, mf),
        nobs = nrow(mf),
        ncat = split_by_name(ncat, items),
        level = split_by_name(level, items),
        saturated = saturated_model,
        dimnames = dims)
}

