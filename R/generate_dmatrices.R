#' generate_dmatrices
#'
#' @param files ...
#' @param infofile ...
#'
#' @return ...
generate_dmatrices <- function(files, infofile) {

  info <- utils::read.table(infofile, header = TRUE)

  lfc_data <- NULL
  pval_data <- NULL
  for (i in 1:length(files)){

    data <- utils::read.table(files[i], header = TRUE)

    col_index_lfc <- grep("lfc", colnames(data))
    col_index_pos <- grep("pos.pval", colnames(data))
    col_index_neg <- grep("neg.pval", colnames(data))

    infofile_index <- grep(basename(files[i]), info$filename)
    selection_order <- info$selection[infofile_index]

    colnames(data)[2:ncol(data)] <- paste0(colnames(data)[2:ncol(data)], infofile_index)

    if (selection_order == "pos") {

      if (is.null(lfc_data) & is.null(pval_data)) {
        lfc_data <- data[ , c(1,col_index_lfc)]
        pval_data <- data[ , c(1,col_index_pos)]
      } else {
        lfc_data <- merge(lfc_data, data[ , c(1,col_index_lfc)], by = "Gene", )
        pval_data <- merge(pval_data, data[ , c(1,col_index_pos)], by = "Gene")
      }

    } else if (selection_order == "neg") {

      if (is.null(lfc_data) & is.null(pval_data)) {
        lfc_data <- data[ , c(1,col_index_lfc)]
        pval_data <- data[ , c(1,col_index_neg)]
      } else {
        lfc_data <- merge(lfc_data, data[ , c(1,col_index_lfc)], by = "Gene")
        pval_data <- merge(pval_data, data[ , c(1,col_index_neg)], by = "Gene")
      }

    } else {
      stop('Cannot identify selection order of screen. Check your single-screen
           infofile: Column `selection` must contain "pos", or "neg"!')
    }
  }
  mlist <- list("lfc_data" = lfc_data, "pval_data" = pval_data)
  return(mlist)
}
