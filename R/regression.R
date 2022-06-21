##' @export
regression <- function(object, ...) UseMethod("regression")

##' @export
regression.catlvm <- function(
   object, data = parent.frame(), regression = NULL,
   method = "hybrid", smoothing = TRUE, control = catlvm.control(), ...
) {

}

library(lava)
g <- lvm(eta1 ~ x1+x2)
regression(g) <- c(y1,y2,y3) ~ eta1
latent(g) <- ~eta1
endogenous(g)
exogenous(g)
