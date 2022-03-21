plot.lvFormula <- function(fmla, font = "Helvetica") {
   lv <- sapply(fmla, function(f) {
      if (length(f[[2]]) == 1) deparse(f[[2]])
      else deparse(f[[2]][[1]])
   })
   item <- lapply(fmla, function(f) setdiff(all.vars(f[[3]]), lv))
   path <- sapply(fmla, function(f) {
      if (length(f[[2]]) == 1) lv <- deparse(f[[2]])
      else lv <- deparse(f[[2]][[1]])
      paste(lv, "-> {", paste(all.vars(f[[3]]), collapse = ", "), "}")
   })
   cl <- sapply(fmla, function(f)
      if (length(f[[2]]) == 3) deparse(f[[2]][[3]])
      else NA)

   txt_def <- paste0(
      "node [shape = box]\n",
      paste(unlist(item), collapse = ", "), "\n\n",
      "node [shape = oval]\n",
      paste(lv, collapse = ", ")
   )
   txt_path <- paste0(path, collapse = "\n")

   if (sum(!is.na(cl)) > 0) {
      txt_grp <- paste0(sapply(which(!is.na(cl)), function(x)
         paste0("subgraph cluster_", cl[x],
                "{ \nlabel = 'in ", cl[x],
                "'\nlabelloc = 'b' \n",
                paste(lv[x], collapse = ", "), "\n",
                paste(item[[x]], collapse = ", "), "\n}"
                )), collapse = "\n\n")
   }

   text <- paste0(
      "digraph { \n", "node[fontname = '", font, "']",
      txt_def, "\n\n", txt_path, "\n\n",
      if (sum(!is.na(cl)) > 0) txt_grp, "\n}"
   )
   # cat(text)
   grViz(text)
}
