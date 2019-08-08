#' BLCAsubtyping - Predicting bladder cancer molecular subtypes
#' 
#' This package provides the 'classify' function, that can be used to predict class for bladder 
#' cancer samples using several differente classification systems, from the sample's expression profile.
#' 
#' @seealso \link[BLCAsubtyping]{classify}
#' @author Aurelien de Reynies, Aurelie Kamoun
#' @note This is a contribution from the Tumor Identity Cards (CIT) program founded by the 
#' 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. 
#' For any question please contact \url{CITR@ligue-cancer.net}
"_PACKAGE"

#' CIT BLCA cohort data
#' 
#' Gene expression data derived from the CIT cohort of bladder cancer.
#' 
#' @format A data.frame that contains gene expression data derived from microarray, 
#'     with probes in the rows and samples in the columns.
#' @references \url{https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-1803/}
"cit"