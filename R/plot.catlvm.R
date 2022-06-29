#' @export
plot.catlvm <- function(x, abbreviation = FALSE, dir = "TD",
                        font = "Helvetica", ...) {
   if (abbreviation) {
      manifest = lapply(x$model$vars$manifest, function(x)
         paste0("'", x[1], " ~ ", x[length(x)], "'"))
   } else {
      manifest = lapply(x$model$vars$manifest, function(x)
         paste0("'", x, "'"))
   }
   latent = lapply(x$model$vars$latent, function(x)
      paste0("'", x, "'"))

   var_def <- paste0(
      "node [shape = box]\n",
      paste(unlist(manifest), collapse = ", "),
      "\n\n node [shape = oval]\n",
      paste(c(x$model$label), collapse = ", ")
   )
   msr_path = paste(sapply(names(manifest), function(x)
      paste(x, "-> {", paste(manifest[[x]], collapse = ", "), "}")),
      collapse = "\n")
   str_path <- paste(sapply(names(latent), function(x)
      paste(x, "-> {", paste(latent[[x]], collapse = ", "), "}")),
      collapse = "\n")

   text <- paste0(
      "digraph { \n  rankdir = '", dir, "';",
      "node[fontname = '", font, "']\n\n",
      var_def, "\n\n", msr_path, "\n",
      str_path, "\n\n}"
   )

   DiagrammeR::grViz(text, ...)
}


