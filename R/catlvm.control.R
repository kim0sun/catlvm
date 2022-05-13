##' @export
catlvm.control <- function(
   em.iterlim = 1000, em.tol = 1e-5,
   nlm.iterlim = 500, nlm.tol = 1e-6,
   verbose = TRUE, per.iter = 1000,
   init.param = NULL, contraints = NULL, ...
) {
   ctrl <- list(
      em.iterlim = em.iterlim, em.tol = em.tol,
      nlm.iterlim = nlm.iterlim, nlm.tol = nlm.tol,
      verbose = verbose, per.iter = per.iter,
      init.param = init.param, contraints = contraints
   )

   class(ctrl) <- "catlvm.control"
   ctrl
}
