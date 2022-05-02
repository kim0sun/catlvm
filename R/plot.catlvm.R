plot.catlvm <- function(x, abbreviation = FALSE, font = "Helvetica", ...) {
   if (abbreviation) {
      manifest = lapply(x$struct$vars$manifest, function(x)
         paste0("'", x[1], " ~ ", x[length(x)], "'"))
   } else {
      manifest = x$struct$vars$manifest
   }
   latent = x$struct$vars$latent

   var_def <- paste0(
      "node [shape = box]\n",
      paste(unlist(manifest), collapse = ", "),
      "\n\n node [shape = oval]\n",
      paste(c(x$struct$label), collapse = ", ")
   )
   msr_path = paste(sapply(names(manifest), function(x)
      paste(x, "-> {", paste(manifest[[x]], collapse = ", "), "}")),
      collapse = "\n")
   str_path <- paste(sapply(names(latent), function(x)
      paste(x, "-> {", paste(latent[[x]], collapse = ", "), "}")),
      collapse = "\n")

   text <- paste0(
      "digraph { \n", "node[fontname = '", font, "']\n\n",
      var_def, "\n\n", msr_path, "\n",
      str_path, "\n\n}"
   )
   DiagrammeR::grViz(text, ...)
   # return(text)
}


