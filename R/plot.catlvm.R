plot.catlvm <- function(x, abbreviation = FALSE, font = "Helvetica", ...) {
   if (abbreviation) {
      manifest = lapply(x$measureModel$manifest, function(x)
         paste0("'", x[1], " ~ ", x[length(x)], "'"))
   } else {
      manifest = x$measureModel$manifest
   }
   var_def <- paste0(
      "node [shape = box]\n",
      paste(unlist(manifest), collapse = ", "),
      "\n\n node [shape = oval]\n",
      paste(unique(c(x$measureModel$label, x$latentStruct$parent)), collapse = ", ")
   )
   msr_path = paste(sapply(names(manifest), function(x)
      paste(x, "-> {", paste(manifest[[x]], collapse = ", "), "}")),
      collapse = "\n")
   str_path <- paste0(x$latentStruct$formula, collapse = "\n")

   text <- paste0(
      "digraph { \n", "node[fontname = '", font, "']\n\n",
      var_def, "\n\n", msr_path, "\n",
      str_path, "\n}"
   )
   DiagrammeR::grViz(text, ...)
   # return(text)
}
