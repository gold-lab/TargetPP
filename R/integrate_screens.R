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
#' @param RRAexact logical. When `TRUE` the RRA algorithm uses their "exact"
#'   method to correct the p-value. Otherwise Bonferroni correction is applied.
#'
#' @return data.frame
#' @export
integrate_screens <- function(datapath,
                              infofile,
                              inputformat = "CASPR",
                              plot = TRUE,
                              RRAexact = TRUE) {

  files <- list.files(datapath, full.names = TRUE)

  # RRA ===============
  gene_lists <- generate_glist(files, infofile)
  RRA <- RobustRankAggreg::aggregateRanks(glist = gene_lists, method = "RRA", exact = RRAexact)
  rownames(RRA) <- NULL

  # RRA pvalue score correction
  if (RRAexact){
    if (min(RRA$Score) < 0){
      warning('The exact p-value correction of the RRA produced negative values! Values are set to the next lowest value...\n
              ("RRA exact p-value correction" can be numerically unstable for small number of screens)' )
      RRA$Score[RRA$Score < 0] <- min(RRA$Score[RRA$Score > 0])
    }
  } else {
    RRA$Score <- RRA$Score * length(gene_lists)
    RRA$Score[RRA$Score > 1] <- 1
  }

  # EBM ===============
  EMP_data <- generate_dmatrices(files, infofile)

  EMP_data_genenames <- EMP_data[["pval_data"]][1]
  EMP_data_pval_cols <- data.matrix(EMP_data[["pval_data"]][-1])
  # data frame with lfc values has to be transposed to be used as the EBM's data
  # matrix. genes have to be the "m variables in rows" and the screens are the
  # "n samples in columns", to assure correct calculation of the data
  # distribution:
  EMP_data_lfc_cols <- t(EMP_data[["lfc_data"]][-1])

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

  # combining p values using the harmonic mean method, but without automatically
  # controlling the FWER in order to subsequently correct for FDR
  data["HMP"] <- apply(data[2:3], MARGIN = 1, FUN = harmonicmeanp::p.hmp, multilevel = FALSE, L = 2)
  data["HMP_fdr"] <- stats::p.adjust(data$HMP, method = "fdr")

  utils::write.csv(data, "valuetable.csv", row.names = F)

  ranking <- data[c(1,5)]
  names(ranking) <- c("Name", "Score")
  ranking <- ranking[order(ranking$Score),]
  utils::write.csv(ranking, "result.csv", row.names = F)

}
