#' @export
summary.catlvm = function(object, type = c("data", "model", "posterior", "parameter"), ...) {
   cat("Summary of CATegorical Latent Variable Model\n")

   if ("data" %in% type) {
      cat("Summary of manifest items\n")
   }
   if ("model" %in% type) {
      cat("Summary of manifest items\n")
   }
   if ("posterior" %in% type) {

   }
   if ("parameter" %in% type) {

   }
}
