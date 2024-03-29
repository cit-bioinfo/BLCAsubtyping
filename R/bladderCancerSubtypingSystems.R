#' classify - Classify a BLCA sample from gene expression data
#' 
#' This function assigns molecular subtypes to bladder samples using expression data, 
#' according to several classification systems.
#' 
#' @param expMat Matrix of normalized expression data. 
#' Samples are in column and probesets are in row. Rownames and colnames are required.
#' @param gpl Dataframe with matching information between expMat rownames and HGNC gene symbols. 
#' The dataframe must have rownames corresponding to ExpMat rownames,
#' and at least one column with associated gene symbols
#' @param symbol Character string referring to the column name of gpl containing gene symbols
#' @param classification.systems A character vector with names of the classification systems to consider
#'  among Baylor, UNC, CIT, Lund, MDA and TCGA.
#'  
#' @return A dataframe of samples annotated according to each classification system.
#' 
#' @author Aurelien de Reynies
#' @keywords methods
#' @note This is a contribution from the Tumor Identity Cards (CIT) program founded by the 
#' 'Ligue Nationale Contre le Cancer' (France): \url{http://cit.ligue-cancer.net}. 
#' For any question please contact \url{CITR@ligue-cancer.net}
#' 
#' @examples
#' data(cit)
#' 
#' #-- Classify using all classification types
#' res <- classify(cit)
#' head(res)
#' 
#' #-- Get TCGA classification
#' tcga_res <- res[, c("ID", "TCGA.subtype")]
#' head(tcga_res)
#' 
#' @import cluster pamr survival
#' @importFrom grDevices graphics.off gray pdf rainbow
#' @importFrom graphics legend par plot
#' @importFrom stats as.dist cor cutree dist hclust mad median pchisq quantile sd setNames t.test var
#' @importFrom utils data
#' @export
#' 
classify <- function(expMat, gpl = NULL, symbol = "Gene.Symbol",
                     classification.systems = c("Baylor", "UNC", "CIT", "Lund", "MDA", "TCGA")
                     ){
  
  if(is.null(gpl)) {
    gpl <- data.frame("Probe.ID" = rownames(expMat), 
                      "Gene.Symbol" = rownames(expMat), 
                      stringsAsFactors = F, row.names = rownames(expMat))
  } else {
    gpl[, symbol] <- as.character(gpl[, symbol])
    if(any(is.na(gpl[, symbol]))) gpl <- gpl[-is.na(gpl[, symbol]), ]
  }
  if(!("matrix" %in% class(expMat))) {
    stop("`expMat` must be a matrix.")
   }
  
  res.class <- data.frame(ID = colnames(expMat), stringsAsFactors = F)
  
  if("Baylor" %in% classification.systems){
    cat("predicting Baylor subtypes...")
    res.class$Baylor.subtype <- baylor.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("...DONE \n")
  }
  
  if("UNC" %in% classification.systems){
    cat("predicting UNC subtypes...")
    res.class$UNC.subtype <- chapelHill.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("\n...DONE \n")
  }
  
  if("CIT" %in% classification.systems){
      cat("predicting CIT subtypes...")
      res.class$CIT.subtype <- cit.classify(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
      cat("...DONE \n")
  }
  
  if("Lund" %in% classification.systems){
      cat("predicting Lund subtypes...")
      res.class$Lund.subtype <- lund.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
      cat("...DONE \n")
  }
  
  if("MDA" %in% classification.systems){
    cat("predicting MDA subtypes...")
    res.class$MDA.subtype <- MDA.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("...DONE \n")
  }
  
  if("TCGA" %in% classification.systems){
    cat("predicting TCGA subtypes...")
    res.class$TCGA.subtype <- TCGA.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("...DONE \n")
  }

  return(res.class)
}


baylor.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol"){

  # if(is.null(Gpl)) { # if Gpl is null we expect rownames of Exp to be 'HUGO Gene Symbols'
  #   Gpl <- data.frame("Probe.ID" = rownames(Exp), "Symbol" = rownames(Exp), stringsAsFactors = F, row.names = rownames(Exp))
  # }
  
  G <- intersect(rownames(Exp), rownames(Gpl))
  Exp <- Exp[G, ]
  Gpl <- Gpl[G, ]

  geneID <- rownames(Gpl)[which(Gpl[, Symbol] %in% as.character(baylor.genes[, "Symbol"]))]
  geneExp <- Exp[geneID, ]
  
  hc <- baylor.hcfun(baylor.dfuncc(t(as.matrix(baylor.standardize(geneExp)))))
  hc.clusters <- cutree(hc, k = 2)
  
  tstat <- apply(geneExp[, names(hc.clusters)], 1, function(z) t.test(z~hc.clusters)$statistic)
  symbol2type <- baylor.genes[, 3]
  names(symbol2type) <- baylor.genes[, 1]
  
  if(diff(sapply(split(tstat, symbol2type[Gpl[geneID, Symbol]]),median)) < 0){ 
    CNAM <- c("Basal", "Differentiated")
  }else{
    CNAM <- c("Differentiated", "Basal")}
  
  #print(split(tstat, symbol2type[Gpl[geneID,Symbol]]))

  Baylor.subtype <- CNAM[hc.clusters]
  names(Baylor.subtype) <- names(hc.clusters)
  
  return(Baylor.subtype)

}

chapelHill.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol", nmin = 10, get_posterior = F){
  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol)
  
  classes <- chapelHill.training$classes
  expMat  <- chapelHill.training$expMat
  genes   <- chapelHill.training$genes
                 
  G <- intersect(genes, rownames(Exp))
  cat(paste(length(G), " of ", nrow(expMat), " genes from the initial predictor are measured in this dataset\n", sep=""))
  
  if(length(G) < nmin) stop("Too few genes to get a prediction.")
  Exp <- Exp[G, ]
  
  trainSet <- list(x = scale(expMat[G, ]), y = classes, geneid = G, genenames = G) #create the training set reduced to common genes
  mytrain  <- pamr.train(trainSet) #train PAM on training dataset

  ChapelHill.subtype <- c("Basal","Luminal")[pamr.predict(mytrain, scale(Exp), 0)] #PAM class predictions
  names(ChapelHill.subtype) <- colnames(Exp)
  
  if(get_posterior) {
    pred.prob  <- pamr.predict(mytrain, scale(Exp), 0, type = "posterior") #PAM prediction probabilities
    colnames(pred.prob) <- c("ChapelHill.posterior.Basal", "ChapelHill.posterior.Luminal")
    ChapelHill.subtype <- data.frame(id = colnames(Exp), ChapelHill.subtype, pred.prob, stringsAsFactors = F)
  }
  
  return(ChapelHill.subtype)
}


cit.classify <- function(Exp, Gpl=NULL, Symbol="Symbol"){
  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol, classification.system = "CIT")
  G <- intersect(rownames(pred[[1]]$mean), rownames(Exp))
  D0  <- pred[[1]]$centroidsdata$data
  CL0 <- pred[[1]]$centroidsdata$samples
  
  d.scale <- D0[G, CL0[, 1]]

  pred. <- cit.centroids(d.scale, CL0[, 2], rowCentering = NULL, dist.meth = "pearson")
  centro.scale <- pred.$centroids$mean
  d.pred <- Exp[G, ]
      
  pred.test <- apply(cor(d.pred, centro.scale, use = "complete.obs"), 1, function(x){colnames(centro.scale)[which.max(x)]})
      
  CIT.subtype <- paste("MC", 1:7, sep = "")[match(pred.test, paste("CIT.CC", c(2, 1, 3:7), sep = ""))]
  
  names(CIT.subtype) <- colnames(Exp)
      
  return(CIT.subtype)
}


lund.predict <- function(Exp, Gpl=NULL, Symbol="Symbol", rowcentering=TRUE, nmin=100, ExplicitClassName = T){
  
  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol, classification.system = "Lund")
  Exp <- t(scale(t(Exp), scale = F))
  Lund.subtype <- ncc.corr(lund.centroids, Exp)
  Lund.subtype[which(Lund.subtype %in% c("GU-Inf1", "GU-Inf2"))] <- "GU-Inf"
  
  return(Lund.subtype)
}


MDA.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol", nmin = 2000){
  
  if(is.null(Gpl)) {
    Gpl <- as.data.frame(cbind("Probe.ID"= rownames(Exp), "Symbol"= rownames(Exp)), stringsAsFactors=F, row.names=rownames(Exp))
    names(Gpl)[2] <- Symbol
  }
  probes <- mda.training$probes
  probes. <- intersect(probes, dimnames(Exp)[[1]])
  d2 <- cit.quantileNormalize(Exp)

  if(length(probes.) > nmin){
    d1 <- mda.training$exp[probes., ]
    d2 <- Exp[probes., ]
  }else{
    G <- intersect(mda.training$gpl[probes, "Symbol"], Gpl[intersect(rownames(Gpl), rownames(Exp)), Symbol])
    probes1 <-  probes[which(mda.training$gpl[probes, "Symbol"] %in% G)]
    sd1 <- apply(mda.training$exp[probes1, ], 1, sd)
    probes1 <- sapply(lapply(split(sd1, mda.training$gpl[probes1, "Symbol"])[-1], which.max), names)
    d1 <- mda.training$exp[probes1, ]
    rownames(d1) <- names(probes1)
    probes2 <- intersect(rownames(Gpl), rownames(Exp))
    probes2 <- probes2[which(Gpl[probes2, Symbol] %in% G)]
    sd2 <- apply(Exp[probes2, ], 1, sd, "na.rm" = T)
    probes2 <- unlist(sapply(lapply(split(sd2, Gpl[probes2, Symbol])[-1], which.max), names))
    d2 <- Exp[probes2, ]
    rownames(d2) <- names(probes2)
  }
      
  G <- intersect(rownames(d1), rownames(d2))
  d1 <- d1[G, ] - apply(d1[G, ], 1, median, "na.rm" = T)
  d2 <- d2[G, ] - apply(d2[G, ], 1, median, "na.rm" = T)
  d <- cbind(d1, d2)
  
  tmp <- as.matrix(dist(t(d)))[colnames(d1), colnames(d2)]
  cl <- mda.training$clin[apply(tmp, 2, function(z) colnames(d1)[which.min(z)]), "Cluster"]
  
  MDA.subtype <- setNames(cl, colnames(d2))

  return(MDA.subtype)

}

TCGA.predict <- function(Exp, Gpl=NULL, Symbol="Symbol"){
  Exp <- cit.probesDataToSymbolData(Exp,Gpl,Symbol,classification.system="TCGA")
  G <- intersect(rownames(tcga.centroids),rownames(Exp))
  
  TCGA.subtype <- colnames(tcga.centroids)[-c(1,2)][apply(cor(Exp[G,], tcga.centroids[G, -c(1,2)], use = "complete.obs"), 1, which.max)]
  names(TCGA.subtype) <- colnames(Exp)
  
  return(TCGA.subtype)

}

