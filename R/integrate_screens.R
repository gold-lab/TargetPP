#' CRISPR multi-screen integration
#'
#' \code{integrate_screens} integrates all CRISPR single-screen results in
#' \code{datapath} and returns a single target-ranking, based on significance
#' and effect size.
#'
#' @param datapath a character string indicating the path to the folder
#'   containing all CRISPR single-screen results.
#' @param infofile a character string indicating the path to the CRISPR
#'   single-screen infofile.
#' @param inputformat a character string indicating the format of the input
#'   data. One of "CASPR" (default), or "MAGeCK".
#' @param plot logical. Should summary plots be generated automatically?
#'
#' @return data.frame
#' @export
integrate_screens <- function(datapath,
                              infofile,
                              inputformat = "CASPR",
                              plot = TRUE) {

  files <- list.files(datapath, full.names = TRUE)

  # RRA ===============
  gene_lists <- generate_glist(files, infofile)
  RRA <- RobustRankAggreg::aggregateRanks(glist = gene_lists, method = "RRA", exact = TRUE)  #TODO: rethink an option to force bonferroni instead of 'exact'

  rownames(RRA) <- NULL

  # EBM ===============
  EMP_data <- generate_dmatrices(files, infofile)

  EMP_data_genenames <- EMP_data[["pval_data"]][1]
  EMP_data_pval_cols <- EMP_data[["pval_data"]][-1]
  EMP_data_lfc_cols <- EMP_data[["lfc_data"]][-1]

  EBM_results <- data.frame(matrix(ncol = 2, nrow = 0))
  colnames(EBM_results) <- c("Name", "EBMscore")


  progbar = utils::txtProgressBar(min = 0, max = nrow(EMP_data_genenames), initial = 0, style = 3)
  for (gene in 1:nrow(EMP_data_genenames)) {
    Name <- EMP_data_genenames[gene,1]
    p_vals_to_aggr <- EMP_data_pval_cols[gene,]
    EBMscore <- EmpiricalBrownsMethod::empiricalBrownsMethod(EMP_data_lfc_cols, p_vals_to_aggr, extra_info = F)

    EBM_results <- rbind(EBM_results, data.frame(Name, EBMscore))

    utils::setTxtProgressBar(progbar,gene)
  }

  EBM <- dplyr::arrange(EBM_results, dplyr::desc(EBMscore))
  colnames(EBM) <- c("Name", "Score")

  rownames(EBM) <- NULL


  # HMP ===============
  data <- merge(RRA, EBM, by="Name")

  names(data) <- c("Name", "RRA_pval", "eBM_pval")

  data["HMP"] <- apply(data[2:3], MARGIN = 1, FUN = harmonicmeanp::p.hmp)  #TODO: investigate whether to set multilevel to FALSE or set L correctly
  data["HMP_fdr"] <- stats::p.adjust(data$HMP, method = "fdr")

  utils::write.csv(data, "valuetable.csv", row.names = F)

  ranking <- data[c(1,5)]
  names(ranking) <- c("Name", "Score")
  ranking <- ranking[order(ranking$Score),]
  utils::write.csv(ranking, "result.csv", row.names = F)

}
