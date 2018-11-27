\name{classify}
\alias{classify}
\title{classify}
\description{Assign molecular subtypes to each bladder samples using expression data, according to several classification systems.}
\usage{classify(expMat, gpl=NULL, symbol="Gene.Symbol",  classification.systems=c("Baylor", "UNC", "CIT", "Lund", "MDA", "TCGA"))}
\arguments{
\item{expMat}{Matrix or dataframe of normalized expression data. Samples are in column and probesets are in row. Rownames and colnames are required.}
\item{gpl}{Dataframe with matching information between expMat rownames and HGNC gene symbols. The dataframe must have rownames corresponding to ExpMat rownames, and at least one column with associated gene symbols}
\item{symbol}{Character string referring to the column name of gpl containing gene symbols}
\item{classification.systems}{A character vector with names of the classification systems to consider among Baylor, UNC, CIT, Lund, MDA and TCGA.}
}
\value{A dataframe of samples annotated according to each classification system.
}
\author{Aurelien de Reynies}
\keyword{methods}
\note{This is a contribution from the Tumor Identity Cards (CIT) program founded by the 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. For any question please contact \url{CITR@ligue-cancer.net}}
