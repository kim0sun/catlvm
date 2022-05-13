proc_data <- function(data, struct) {
   if (is.null(data)) return(NULL)
   f <- paste("~", paste(unlist(struct$vars$manifest), collapse = "+"))
   mf <- model.frame(formula(f), data)
   dims <- dimnames(mf)
   mi <- struct$vars$manifest
   if (any(!unlist(mi) %in% dims[[2]]))
      stop("Following manifest variables not found:\n ",
           paste(unlist(mi)[!unlist(mi) %in% dims[[2]]],
                 collapse = " "))

   data_attr <- lapply(mi, function(x) {
      m <- lapply(mf[x], factor)
      lev <- sapply(m, levels)
      nlev <- sapply(m, nlevels)
      y <- sapply(m, as.numeric)
      rownames(y) <- dims[[1]]
      list(t(y), lev, nlev)
   })

   yf <- do.call(rbind, lapply(data_attr, "[[", 1))
   prl <- prelim.cat(t(yf))
   logpost.cat(prl, em.cat(prl))


   list(nobs = nrow(mf),
        dimnames = dimnames(mf),
        y = unlist(lapply(data_attr, function(x) x[[1]])),
        level = lapply(data_attr, function(x) x[[2]]),
        ncat = lapply(data_attr, function(x) x[[3]]))
}

