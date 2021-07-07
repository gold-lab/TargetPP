#' CRISPR multi-screen integration
#'
#' \code{integrate_screens} integrates all CRISPR single-screen results in
#' \code{datapath} and returns a single target-ranking, based on significance
#' and effect size.
#'
#' @param datapath a character string indicating the path to the CRISPR
#'   single-screen results.
#' @param infofile a character string indicating the path to the CRISPR
#'   single-screen infofile.
#' @param inputformat a character string indicating the format of the input
#'   data. One of "CASPR" (default), or "MAGeCK".
#' @param plot logical. Should summary plots be generated automatically?
#'
#' @return data.frame
#' @export
#'
#' @examples
#' datapath <- ""
#' infofile <- ""
#' integration_result <- integrate_screens(datapath, infofile)
integrate_screens <- function(datapath,
                              infofile,
                              inputformat = "CASPR",
                              plot = TRUE) {

}
