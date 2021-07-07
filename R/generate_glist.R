#' generate_glist
#'
#' @param files ...
#' @param infofile ...
#'
#' @return ...
generate_glist <- function(files, infofile) {
  lfc <- NULL

  info <- utils::read.table(infofile, header = TRUE)

  # generate list of ranklists
  glist <- list()
  for (i in 1:length(files)){

    data <- utils::read.table(files[i], header = TRUE)
    col_index <- grep("lfc", colnames(data))

    genenames_and_values <- data[ , c(1,col_index)]

    infofile_index <- grep(basename(files[i]), info$filename)
    selection_order <- info$selection[infofile_index]

    if (selection_order == "pos") {
      ordered_by_col_values <- dplyr::arrange(genenames_and_values, dplyr::desc(lfc))  #TODO: check if this should be the other way around
    } else if (selection_order == "neg") {
      ordered_by_col_values <- dplyr::arrange(genenames_and_values, lfc)
    } else {
      stop('Cannot identify selection order of screen. Check your single-screen
           infofile: Column `selection` must contain "pos", or "neg"!')
    }
    glist[length(glist)+1] <- list(as.character(ordered_by_col_values[[1]]))
  }
  return(glist)
}
