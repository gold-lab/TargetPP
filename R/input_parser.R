#' Parsing inputs to genelists
#'
#' Parsing input files to the genelist format used by the RobustRankAggregation.
#'
#' @param files vector of character strings representing paths to all CRISPR
#'   single-screen results.
#' @param infotable data.frame containing information about the individual
#'   CRISPR screen results.
#'
#' @return list of genelists
parse_glist <- function(files, infotable) {
  lfc <- NULL

  # generate list of ranklists
  glist <- list()
  for (i in 1:length(files)){

    data <- utils::read.table(files[i], header = TRUE)
    col_index <- grep("lfc", colnames(data))

    genenames_and_values <- data[ , c(1,col_index)]

    infotable_index <- grep(basename(files[i]), infotable$filename)
    selection_order <- infotable$selection[infotable_index]

    if (selection_order == "pos") {
      ordered_by_col_values <- dplyr::arrange(genenames_and_values,
                                              dplyr::desc(lfc))
    } else if (selection_order == "neg") {
      ordered_by_col_values <- dplyr::arrange(genenames_and_values,
                                              lfc)
    } else {
      stop('Cannot identify selection order of screen. Check your single-screen
           infofile: Column `selection` must contain "pos", or "neg"!')
    }
    glist[length(glist)+1] <- list(as.character(ordered_by_col_values[[1]]))
  }
  return(glist)
}



#' Parsing inputs to data matricies
#'
#' Parsing input files to the matrix format used by the empirical Brown's
#' method.
#'
#' @param files vector of character strings representing paths to all CRISPR
#'   single-screen results.
#' @param infotable data.frame containing information about the individual
#'   CRISPR screen results.
#'
#' @return list of two matricies
parse_dmatrices <- function(files, infotable) {
  lfc_data <- NULL
  pval_data <- NULL

  for (i in 1:length(files)){

    data <- utils::read.table(files[i], header = TRUE)

    col_index_lfc <- grep("lfc", colnames(data))
    col_index_pos <- grep("pos.pval", colnames(data))
    col_index_neg <- grep("neg.pval", colnames(data))

    infotable_index <- grep(basename(files[i]), infotable$filename)
    selection_order <- infotable$selection[infotable_index]

    colnames(data)[2:ncol(data)] <- paste0(colnames(data)[2:ncol(data)],
                                           infotable_index)

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
