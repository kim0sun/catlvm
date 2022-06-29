#' @export
confint.catlvm <- function(
   object, parm, level = 0.95,
   method = c("asymp", "logit"),
   out = c("param", "logit")
) {
   if (missing(parm)) parm <- c("pi", "tau", "rho")
   out <- match.arg(out)

   logit <- unlist(object$estimates$logit[parm])
   se <- unlist(object$estimates$se$logit[parm])

   lower <- (1 - level) / 2
   upper <- 1 - lower
   cn <- format.pc(c(lower, upper), 3)

   ci_logit <- logit + se %o% qnorm(c(lower, upper))
   args <- object$args

   switch(out, param = {
      ci_par <- sapply(seq_len(ncol(ci_logit)), function(i) {
         cii <- logit2log(
            ci_logit[, i], args$ncat, args$nroot,
            args$nlink_unique, args$nleaf_unique,
            args$root - 1, args$u - 1, args$v - 1,
            args$nclass_root, args$nclass_leaf,
            args$nclass_u, args$nclass_v)
         unlist(output_param(cii, object$model, args))
      })
      ci <- t(apply(ci_par, 1, sort))
      colnames(ci) <- cn
      ci
   }, logit = {
      ci <- sapply(seq_len(ncol(ci_logit)), function(i) {
         cii <- splitlogit(
            ci_logit[, i], args$ncat, args$nroot,
            args$nlink_unique, args$nleaf_unique,
            args$root, args$u, args$v, args$nclass_root,
            args$nclass_u, args$nclass_v, args$nclass_leaf
         )
         unlist(output_logit(cii, object$model, args))
      })
      colnames(ci) <- cn
      ci
   })
}

format.pc <- function(perc, digits)
   paste(format(100 * perc, trim = TRUE, scientific = FALSE, digits = digits), "%")
