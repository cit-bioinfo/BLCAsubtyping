# Rd
# description >> Assign molecular subtypes to each bladder samples using expression data, according to several classification systems.
# argument
# item >> expMat >> Matrix of expression data
# item >> gpl >> required
# item >> symbol >> required
# item >> classification.systems >> A character vector with names of the classification systems to consider among Baylor, UNC, CIT, Lund, MDA and TCGA.
# value >> A dataframe of samples annotated according to each classification system.
# author >> Aurelien de Reynies
# keyword >> methods
# details >> ...
# seealso >> ...
# references >> ...
# examples >> ...
# end
classify <- function(expMat, gpl = NULL, symbol = "Gene.Symbol",  
                     classification.systems = c("Baylor", "UNC", "CIT", "Lund", "MDA", "TCGA")
                     ){
  
  res.class <- data.frame(ID = colnames(expMat), stringsAsFactors = F)
  
  if("Baylor" %in% classification.systems){
    cat("predicting Baylor subtypes...")
    res.class$Baylor.subtype <- baylor.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("...DONE \n")
  }
  
  if("UNC" %in% classification.systems){
    cat("predicting UNC subtypes...")
    res.class$UNC.subtype <- chapelHill.predict(expMat, Gpl = gpl, Symbol = symbol)[res.class$ID]
    cat("...DONE \n")
  }
  
  if("CIT" %in% classification.systems){
      cat("predicting CIT subtypes...")
      res.class$CIT.subtype <- cit.classify(expMat, Gpl = gpl, Symbol = symbol, classif = "CC")[res.class$ID]
      cat("...DONE \n")
  }
  
  if("Lund" %in% classification.systems){
      cat("predicting Lund subtypes...")
      res.class$Lund.subtype <- lund.predict(expMat, Gpl = gpl, Symbol = symbol, classif = "Lund2017")[res.class$ID]
      cat("...DONE \n")
  }
  
  if("MDA" %in% classification.systems){
    cat("predicting MDA subtypes...")
    res.class$MDA.subtype <- MDA.predict(expMat, Gpl = gpl, Symbol = symbol, distance.method = "euclid")[res.class$ID]
    cat("...DONE \n")
  }
  
  if("TCGA" %in% classification.systems){
    cat("predicting TCGA subtypes...")
    res.class$TCGA.subtype <- TCGA.predict(expMat, Gpl = gpl, Symbol = symbol, classif = "2017")[res.class$ID]
    cat("...DONE \n")
  }

  return(res.class)
}


baylor.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol"){
  
  data(baylor.genes)
      
  if(is.null(Gpl)) { # if Gpl is null we expect rownames of Exp to be 'HUGO Gene Symbols'
    Gpl <- as.data.frame("Probe.ID" = rownames(Exp), "Symbol" = rownames(Exp), stringsAsFactors = F, row.names = rownames(Exp))
  }
  
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
      
  require(pamr)
  data(chapelHill.training)
  
  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol)
  
  classes <- chapelHill.training$classes
  expMat  <- chapelHill.training$expMat
  genes   <- chapelHill.training$genes
                 
  G <- intersect(genes, rownames(Exp))
  print(paste(length(G), " of ", nrow(expMat), " genes from the initial predictor are measured in this dataset", sep=""))
  
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


cit.classify <- function(Exp, Gpl=NULL, Symbol="Symbol", classif = c("BL", "CC")[1]){
  
  if(classif == "BL") data(CIT.BL_predictor) else if(classif == "CC") data(CIT.CC_predictor) else if(classif == "cs98") data(CIT.cs98_predictor)

  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol, classification.system = "CIT")
  G <- intersect(rownames(pred[[1]]$mean), rownames(Exp))
  D0  <- pred[[1]]$centroidsdata$data
  CL0 <- pred[[1]]$centroidsdata$samples
  
  if(classif == "BL"){
      
      pred. <- cit.centroids(D0[G,CL0[,1]], CL0[, 2], median, dist.meth ="pearson")
      p <- cit.distToCentroids(Exp[G,], pred.[[1]], dist.meth="pearson", d.isPretreated = T)
      
      CIT.subtype <- c("Other", "InvBL")[as.numeric(p[[4]] == "Inv_BL") + 1]
  }
  
  if(classif == "CC"){
      d.scale <- D0[G, CL0[, 1]]

      pred. <- cit.centroids(d.scale, CL0[, 2], rowCentering = NA, dist.meth = "pearson")
      centro.scale <- pred.$centroids$mean
      d.pred <- Exp[G, ]
      
      pred.test <- apply(cor(d.pred, centro.scale, use = "complete.obs"), 1, function(x){colnames(centro.scale)[which.max(x)]})
      
      CIT.subtype <- paste("MC", 1:7, sep = "")[match(pred.test, paste("CIT.CC", c(2, 1, 3:7), sep = ""))]
      
  }
  names(CIT.subtype) <- colnames(Exp)
      
  return(CIT.subtype)
}


lund.predict <- function(Exp, Gpl=NULL, Symbol="Symbol", classif = "Lund2017", rowcentering=TRUE, nmin=100, ExplicitClassName = T){
  
  Exp <- cit.probesDataToSymbolData(Exp, Gpl, Symbol, classification.system = "Lund")
  
  if(classif == "Lund2012"){
      
      data(lund.centroids)
      
      G <- intersect(rownames(lund.centroids), rownames(Exp))
      
      if(length(G) < nmin) stop("Too few genes to get a prediction.")
      
      centroids <- lund.centroids[G, ]
      Exp <- Exp[G, ]
      
      if(rowcentering) Exp <- Exp - apply(Exp, 1, median, "na.rm" = T)
      
      tmp <-  as.data.frame(cor(Exp, centroids, use="pair"))
      names(tmp) <- paste("corTo", names(centroids), sep="")
      
      Lund.subtype <-  names(centroids)[apply(tmp, 1, which.max)]
      
      if (ExplicitClassName) {
          Lund.subtype <- c("Uro.A", "Uro.A", "Geno.Unst", "Geno.Unst", "Infilt", "Uro.B", "SCC-l")[match(
          Lund.subtype,c("MS1a", "MS1b", "MS2a1", "MS2a2", "MS2b1", "MS2b2.1", "MS2b2.2"))]
      }
      names(Lund.subtype) <- colnames(Exp)

  }
  
  if(classif == "Lund2017"){
      
      data(lund2017.centroids)
      
      Exp <- t(scale(t(Exp), scale = F))
      Lund.subtype <- ncc.corr(lund.centroids, Exp)
      Lund.subtype[which(Lund.subtype %in% c("GU-Inf1", "GU-Inf2"))] <- "GU-Inf"
      
  }
  
  return(Lund.subtype)
}


MDA.predict <- function(Exp, Gpl = NULL, Symbol = "Symbol", nmin = 2000, distance.method = c("euclid", "pearson")[2]){
  
  if(is.null(Gpl)) {
    Gpl <- as.data.frame(cbind("Probe.ID"= rownames(Exp), "Symbol"= rownames(Exp)), stringsAsFactors=F, row.names=rownames(Exp))
    names(Gpl)[2] <- Symbol
  }

  data(mda.training)

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
  
  if(distance.method == "euclid"){
    tmp <- as.matrix(dist(t(d)))[colnames(d1), colnames(d2)]
    cl <- mda.training$clin[apply(tmp, 2, function(z) colnames(d1)[which.min(z)]), "Cluster"]
  } else if (distance.method == "pearson"){
    tmp <- cor(d1, d2, use = "complete.obs", method = "pearson")
    cl <- mda.training$clin[apply(tmp, 2, function(z) colnames(d1)[which.max(z)]), "Cluster"]
    #clusters <- names(table(mda.training$clin$Cluster))
    #cl <- clusters[unlist(apply(aggregate(tmp, by = list(mda.training$clin[, "Cluster"]), median)[, -1], 2, which.max))]
  }
  
  MDA.subtype <- setNames(cl, colnames(d2))

  return(MDA.subtype)

}



TCGA.predict <- function(Exp, Gpl=NULL, Symbol="Symbol", classif = c("2017", "2013")[1]){
  
  if(classif == "2017") data(tcga2017.centroids) else data(tcga.centroids)
  
  Exp <- cit.probesDataToSymbolData(Exp,Gpl,Symbol,classification.system="TCGA")
  G <- intersect(rownames(tcga.centroids),rownames(Exp))
  
  TCGA.subtype <- colnames(tcga.centroids)[-c(1,2)][apply(cor(Exp[G,], tcga.centroids[G, -c(1,2)], use = "complete.obs"), 1, which.max)]
  names(TCGA.subtype) <- colnames(Exp)
  
  return(TCGA.subtype)

}

