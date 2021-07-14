#' CRISPR multi-screen integration
#'
#' \code{integrate_screens} integrates all CRISPR single-screen results in
#' \code{datapath} and returns a single target-ranking, based on significance
#' and effect size across all input screens.
#'
#' @param datapath a character string indicating the path to the folder
#'   containing all CRISPR single-screen results.
#' @param infofile a character string indicating the path to the CRISPR
#'   single-screen infofile.
#' @param input_format a character string indicating the format of the input
#'   data. Must be one of "CASPR" (default), or "MAGeCK".
#' @param plot logical. Should summary plots be generated automatically?
#' @param RRAexact logical. When `TRUE` (default) the RRA algorithm uses its
#'   "exact" method to correct the p-value. Otherwise Bonferroni correction is
#'   applied.
#' @param extra_info logical. When `FALSE` (default) only the adjusted P-value
#'   is returned. Setting this to `TRUE` adds columns for each intermediary
#'   value to the returned data.frame.
#'
#' @return data.frame
#' @export
integrate_screens <- function(datapath,
                              infofile,
                              input_format = "CASPR",
                              extra_info = FALSE,
                              plot = FALSE,
                              RRAexact = TRUE) {
  HMP_padj <- NULL

  files <- list.files(datapath, full.names = TRUE)
  infotable <- utils::read.table(infofile, header = TRUE)

  # RobustRankAggregation ======================================================
  gene_lists <- parse_glist(files, infotable)
  RRA <- RobustRankAggreg::aggregateRanks(glist = gene_lists,
                                          method = "RRA",
                                          exact = RRAexact)
  rownames(RRA) <- NULL

  # RRA pvalue score correction
  if (RRAexact){
    if (min(RRA$Score) < 0){
      warning('The exact p-value correction of the RRA produced negative values!
              Values are set to the next lowest value...\n
              ("RRA exact p-value correction" can be numerically unstable for
              small number of screens)' )
      RRA$Score[RRA$Score < 0] <- min(RRA$Score[RRA$Score > 0])
    }
  } else {
    RRA$Score <- RRA$Score * length(gene_lists)
    RRA$Score[RRA$Score > 1] <- 1
  }
  colnames(RRA) <- c("Gene", "RRA_pval")

  # Empirical Brown's Method ===================================================
  EMP_data <- parse_dmatrices(files, infotable)

  EMP_data_genenames <- EMP_data[["pval_data"]][1]
  EMP_data_pval_cols <- data.matrix(EMP_data[["pval_data"]][-1])
  # data frame with lfc values has to be transposed to be used as the EBM's data
  # matrix. genes have to be the "m variables in rows" and the screens are the
  # "n samples in columns", to assure correct calculation of the data
  # distribution:
  EMP_data_lfc_cols <- t(EMP_data[["lfc_data"]][-1])

  EBM <- data.frame(matrix(ncol = 2, nrow = 0))

  progbar = utils::txtProgressBar(min = 0,
                                  max = nrow(EMP_data_genenames),
                                  initial = 0,
                                  style = 3)
  for (gene in 1:nrow(EMP_data_genenames)) {
    genename <- EMP_data_genenames[gene,1]
    p_vals_to_aggr <- EMP_data_pval_cols[gene,]
    EBMscore <- EmpiricalBrownsMethod::empiricalBrownsMethod(EMP_data_lfc_cols,
                                                             p_vals_to_aggr,
                                                             extra_info = FALSE)

    EBM <- rbind(EBM, data.frame(genename, EBMscore))

    utils::setTxtProgressBar(progbar,gene)
  }
  colnames(EBM) <- c("Gene", "EBM_pval")
  rownames(EBM) <- NULL

  # HMP ===============
  TPP_result <- merge(RRA, EBM, by="Gene")

  # combining p values using the harmonic mean method, but without automatically
  # controlling the FWER, in order to subsequently correct for FDR:
  TPP_result["HMP_pval"] <- apply(TPP_result[2:3],
                       MARGIN = 1,
                       FUN = harmonicmeanp::p.hmp,
                       multilevel = FALSE,
                       L = 2)
  TPP_result["HMP_padj"] <- stats::p.adjust(TPP_result$HMP_pval, method = "fdr")

  TPP_result <- dplyr::arrange(TPP_result, HMP_padj)

  if (!extra_info) { TPP_result <- TPP_result[c(1,5)] }
  return(TPP_result)

}
