plotSurv <- function(Clin, Event="OS", Time="OSdelay", File=NULL,classification.system=c("Baylor","ChapelHill","CIT","Lund","MDA","TCGA")[1],posLegend="bottom",main.prefix=""){
     require(survival) 
     
     if(classification.system=="Baylor"){
          Strat="Baylor.subtype"
          StratColors=c("Basal"="red3","Differentiated"="cyan")     
     }
     
     if(classification.system=="ChapelHill"){
          Strat="ChapelHill.subtype"
          StratColors=c("Basal"="red3","Luminal"="cyan")
     }
     
     if(classification.system=="CIT"){
         Strat="CIT.subtype"
         StratColors=c("basal-like"="red3","other"="cyan")     
     }
     
     if(classification.system=="Lund"){
         Strat="Lund.subtype2" 
         StratColors=c("MS1a"="darkseagreen","MS1b"="green3","MS2a1"="lightskyblue","MS2a2"="royalblue3","MS2b1"="pink","MS2b2.1"="saddlebrown","MS2b2.2"="red2","Uro.A"="darkseagreen",  "Geno.Unst"="lightskyblue", "Infilt"="pink", "Uro.B"="saddlebrown", "SCC-l"="red2")
         StratColors <- StratColors[unique(as.character(Clin[,Strat]))]
     }
     
     if(classification.system=="MDA"){
          Strat="MDA.subtype"
          StratColors=c("basal"="red3","luminal"="green","p53-like"="cyan")
     }
     
     if(classification.system=="TCGA"){
          Strat="TCGA.subtype"
          StratColors=c("class 1"="green","class 2"="cyan","class 3"="red3","class 4"="orange")
     }

      if(!is.null(File)) pdf(file=File,height=8.5,width=11.5)
      chisq <-pval<-NA
      try(tmp <- survdiff(Surv(Clin[,Time],Clin[,Event]) ~ Clin[,Strat]))
      try(chisq <- tmp$chisq )
      try(df    <- (sum(1 * (tmp$exp > 0))) - 1 )
      try(pval <- pchisq(chisq,df,lower.tail=F) )

      par(oma=c(2,3,0,0),mar=c(5,5,3,3))
      try(fitOS<-survfit(Surv(Clin[,Time],Clin[,Event])~ Clin[,Strat]))
      try(plot(fitOS,cex=1.5,cex.lab=1.5,cex.axis=1.5,ylab="Probability",xlab="Months",col=StratColors,main=paste(main.prefix,"Chi2= ",round(chisq,digits=2),", df=",df,", log-rank p-value=",format(pval,digits=3),sep=""),lwd=3))
      try(legend(posLegend,legend=names(StratColors),col=StratColors,lty=1,bty="n",cex=2,lwd=3))
      if(!is.null(File))graphics.off()
}




cit.probesDataToSymbolData <- function(expMat,gpl=NULL,symbol="Symbol",opt.fun=NULL,more=FALSE,classification.system=c("CIT","Lund","TCGA")[1]){
         if(classification.system=="CIT")opt.fun <- NULL
         if(classification.system=="Lund")opt.fun <- cit.robustCV
         if(classification.system=="TCGA")opt.fun <- mad
         

         if(is.null(gpl)) return(expMat) # if gpl is null then the rownames of the expression matrix are expected to be already HUGO Gene Symbols 

         # these lines are to the test the situation of a non null gpl but with an expression matrix which is already aggregated by HUGO Gene Symbols
         rn <- intersect(rownames(expMat),rownames(gpl))
         rs <- intersect(rownames(expMat),gpl[,symbol])
         if(length(rs) > length(rn)) return(expMat)

         if(is.null(opt.fun)) return(cit.dfAggregate(expMat[rn,],gpl[rn,symbol]))
         
         val <- apply(expMat[rn,],1,opt.fun)
         names(val) <- rn

         selectedProbes <- unlist(sapply(split(val,gpl[rn,symbol]),function(z)names(z)[which.max(z)]))

         expSymbol <- expMat[selectedProbes,]
         rownames(expSymbol) <- gpl[selectedProbes,symbol]

         if(more){
              list(expSymbol=expSymbol,probesToSymbol=gpl[selectedProbes,symbol,drop=F])
         }else{
              expSymbol
         }

}


cit.dfAggregate <- function (data, partition, MARGIN = 1, fAggreg = mean.na)    {

    if (MARGIN == 1 & (identical(fAggreg, mean.na) | identical(fAggreg, 
        mean) | identical(fAggreg, sum) | identical(fAggreg, 
        sum.na)) & !any(is.na(colMeans(data)))) {
        res <- rowsum(data, partition)
        if (identical(fAggreg, mean.na) | identical(fAggreg, 
            mean)) {
            len <- sapply(split(1:nrow(data), partition), length)
            res <- res/len[rownames(res)]
        }
        return(res)
    }
    cMARGIN <- setdiff(c(1, 2), MARGIN)
    n <- length(partition)
    N <- dim(data)[MARGIN]
    p <- dim(data)[cMARGIN]
    if (n != N) 
        stop("Error - function cit.dfAggregate : size of partition doesn't correspond to data dimension")
    l <- split(1:N, partition)
    d <- data
    if (MARGIN == 2) 
        d <- t(data)
    d <- matrix(sapply(l, function(i) if (length(i) == 1) {
        unlist(d[i, ])
    }
    else {
        apply(d[i, ], 2, fAggreg)
    }), ncol = p, byrow = TRUE)
    d <- as.data.frame(d)
    rownames(d) <- names(l)
    names(d) <- dimnames(data)[[cMARGIN]]
    if (MARGIN == 2) 
        d <- as.data.frame(t(d))
    d
}

mean.na <- function(z){mean(z,na.rm=TRUE)}
median.na <- function(z){median(z,na.rm=TRUE)}
var.na <- function(z){var(z,na.rm=TRUE)}
sum.na <- function(z){sum(z,na.rm=TRUE)}

cit.rainbow  <- function (n, default2Col = c("black", "red"), default6Col = c("red",
                          "blue", "green", "orange", "purple", "gray"), colorMode = c("rainbow",
                          "gray", "grey")) {
    colorMode = match.arg(colorMode)
    if (colorMode == "grey")
        colorMode = "gray"
    if (n == 1)
        return(default2Col[1])
    if (n == 2)
        return(default2Col)
    if (n > 6 & colorMode == "rainbow")
        return(rainbow(n))
    if (n < 7 & n > 2 & colorMode == "rainbow")
        return(default6Col[1:n])
    if (n > 2 & colorMode == "gray")
        return(gray(seq(1, 0, length = n)))
}

cit.nomValToColor  <- function (val, col = NULL, colna = "grey", colpos = "red",
                                colneg = "gray", colfun = cit.rainbow){
    tval <- table(val)
    vals <- names(tval)
    n <- length(tval)
    if (is.null(col)) {
        if (n == 2) {
            col <- c(colneg, colpos)
            neg <- intersect(vals, c("WT", "Wt", "wt", "mss",
                "MSS", "-", "0", "N", "n", "NEG", "Neg", "No",
                "NO", "NON", "Non", "non"))
            if (length(neg) > 0) {
                names(col) <- c(neg, setdiff(vals, neg))
            }
        }
        else {
            col <- colfun(n)
        }
    }
    if (!is.null(names(col))) {
        valToCol <- col[as.character(val)]
    }
    else {
        valToCol <- col[as.numeric(as.factor(val))]
    }
    valToCol[which(is.na(valToCol))] <- colna
    list(valToCol = valToCol, col = col)
}




baylor.standardize = function(mymatrix){
  emc.rm = apply(mymatrix,1,mean,na.rm=TRUE)
  emc.rsd = apply(mymatrix,1,sd,na.rm=TRUE)
  emcmicnorm = sweep(mymatrix,1,emc.rm,FUN="-")
  emcmicnorm = sweep(emcmicnorm,1,emc.rsd,FUN="/")
  emcmicnorm
}

baylor.dfuncc = function(d){
  td = t(d)
  ds= as.dist(1-cor(td))
  attr(ds, "Size") <- ncol(td)
  attr(ds, "Labels") <- colnames(td)
  attr(ds, "Diag") <- FALSE
  attr(ds, "Upper") <- FALSE
  attr(ds, "method") <- "1-cor(d)"
  return(ds)
}

baylor.hcfun = function(d){hclust(d, method="complete")}



cit.robustCV <- function(z){
    mad(z,"na.rm"=T)/median(z,"na.rm"=T)
}


cit.quantileNormalize <- function (d, baseline = NULL) {
    nr <- nrow(d)
    if (is.null(baseline)) 
        baseline <- quantile(as.vector(d), probs = seq(0.5, nr - 
            0.5, by = 1)/nr, na.rm = TRUE)
    for (i in 1:ncol(d)) {
        ok <- which(!is.na(d[, i]))
        nOk <- length(ok)
        if (nOk == nr) {
            d[order(d[, i]), i] <- baseline
        }
        else {
            if (nOk > 1) {
                d[ok[order(d[ok, i])], i] <- quantile(baseline, 
                  probs = seq(0.5, nOk - 0.5, by = 1)/nOk)
            }
            else {
                stop("Some columns have to few non NA values : these columns should be removed.")
            }
        }
    }
    return(d)
}



cit.centroids <- function (d, classes, rowCentering = c(NULL, mean, median)[[3]], 
                           rowClassesForAggregation = NULL, rowClassesToKeep = NULL, 
                           dist.meth = c("spearman", "euclidian", "maximum", "manhattan", 
                                         "canberra", "binary", "minkowski", "pearson", "dlda", 
                                         "dqda", "cosine"), maxDist = 0.5, ...) 
{
  n <- length(classes)
  N <- dim(d)[2]
  if (n != N) 
    stop("Error - function cit.centroids: size of classes doesn't correspond to number of columns in data")
  dist.meth <- match.arg(dist.meth)
  if (!is.null(rowClassesForAggregation)) {
    d <- cit.dfAggregate(d, rowClassesForAggregation)
    if (!is.null(rowClassesToKeep)) {
      rowClassesToKeep <- intersect(rowClassesToKeep, rownames(d))
      if (length(rowClassesToKeep) < 2) 
        stop("Error - function cit.centroids  : less than 2 rowClassesToKeep are in d")
      d <- d[rowClassesToKeep, ]
    }
  }
  if (!is.null(rowCentering)) 
    d <- sweep(d, 1, apply(d, 1, rowCentering, na.rm = T))
  w <- which(!is.na(classes))
  if (length(w) == 0) 
    stop("Error - function cit.centroids: classes are undefined (NA)")
  nc <- table(classes)
  m <- cit.dfAggregate(d[, w], classes[w], MARGIN = 2, fAggreg = mean.na)
  v <- cit.dfAggregate(d[, w], classes[w], MARGIN = 2, fAggreg = var.na)
  va <- apply(v, 1, function(vg) sum((nc - 1) * vg)/(length(w) - 
                                                       length(nc)))
  vg <- apply(d[, w], 1, var.na)
  L <- list(mean = m, var = v, `aggregated var` = va, `global var` = vg, 
            rowCentering = rowCentering, rowClassesForAggregation = rowClassesForAggregation, 
            rowClassesToKeep = rowClassesToKeep, centroidsdata = list(samples = data.frame(samplename = colnames(d[, 
                                                                                                                   w]), class = classes[w], stringsAsFactors = F), data = d[, 
                                                                                                                                                                            w]))
  list(centroids = L, distToCentroids = cit.distToCentroids(d, 
                                                            L, dist.meth = dist.meth, maxDist = maxDist, d.isPretreated = TRUE, 
                                                            ...))
}

cit.distToCentroids <- function (d, centroids, dist.meth = c("spearman", "euclidian", 
                                                             "maximum", "manhattan", "canberra", "binary", "minkowski", 
                                                             "pearson", "dlda", "dqda", "cosine"), maxDist = 0.5, d.isPretreated = FALSE, 
                                 sdifftop = NULL, sdisttocent = NULL, verbose = FALSE) 
{
  dist.meth <- match.arg(dist.meth)
  srcs <- colnames(d)
  if (!d.isPretreated) {
    if (!is.null(centroids$rowClassesForAggregation)) {
      d <- cit.dfAggregate(d, centroids$rowClassesForAggregation)
      if (!is.null(centroids$rowClassesToKeep)) {
        centroids$rowClassesToKeep <- intersect(centroids$rowClassesToKeep, 
                                                rownames(d))
        if (length(centroids$rowClassesToKeep) < 2) 
          stop("Error - function cit.distToCentroids: less than 2 rowClassesToKeep are in d")
        d <- d[centroids$rowClassesToKeep, ]
      }
    }
    if (!is.null(centroids$rowCentering)) 
      d <- sweep(d, 1, apply(d, 1, centroids$rowCentering, na.rm = T))
  }
  if (!is.null(centroids$centroidsdata)) {
    dc <- centroids$centroidsdata$data
    colnames(dc) <- paste("centroid.", colnames(dc), sep = "")
    d <- cbind(d, dc)
  }
  N <- ncol(d)
  n <- ncol(centroids$mean)
  if (dist.meth %in% c("dlda", "dqda")) {
    sumlogvar <- apply(log(centroids$var), 2, sum.na)
    if (dist.meth == "dlda") {
      tmp <- apply(d, 2, function(z) apply((z - centroids$mean)^2/centroids$"aggregated var", 
                                           2, sum.na))
      tdist <- as.data.frame(t(tmp))
    }
    if (dist.meth == "dqda") {
      tmp <- apply(d, 2, function(z) apply((z - centroids$mean)^2/centroids$var, 
                                           2, sum.na) + sumlogvar)
      tdist <- as.data.frame(t(tmp))
    }
  }
  else {
    d2 <- t(cbind(d, centroids$mean))
    tdist <- as.matrix(cit.dist(d2, meth = dist.meth, diag = TRUE))
    tdist <- as.data.frame(tdist[1:N, (N + 1):(N + n)])
  }
  rownames(tdist) <- names(d)
  names(tdist) <- names(centroids$mean)
  pred <- apply(tdist, 1, function(z) names(centroids$mean)[which.min(z)])
  mind <- apply(tdist, 1, min)
  difftofirst <- function(x) {
    m <- min(x)
    x - m
  }
  difftop <- apply(tdist, 1, function(x) {
    diff(sort(x)[1:2])
  })
  if (is.null(sdifftop)) {
    if (length(grep("centroid.", names(difftop), value = T)) > 
        0) 
      sdifftop <- round(quantile(difftop[grep("centroid.", 
                                              names(difftop), value = T)], 0.01, na.rm = T), 
                        3)
    else sdifftop <- round(quantile(difftop, 0.01, na.rm = T), 
                           3)
  }
  predf <- sapply(1:nrow(tdist), function(i) {
    if (difftop[i] <= sdifftop) {
      w <- which(difftofirst(tdist[i, ]) <= sdifftop)
      paste(sort(colnames(tdist[i, ])[w]), collapse = "")
    }
    else {
      colnames(tdist[i, ])[which.min(tdist[i, ])]
    }
  })
  names(predf) <- names(pred)
  if (length(grep("centroid.", names(difftop), value = T)) > 
      0) {
    sam <- grep("centroid.", names(difftop), value = T)
    wsure <- sam[centroids$centroidsdata$samples[match(sub("centroid.", 
                                                           "", sam), centroids$centroidsdata$samples$samplename), 
                                                 "class"] == predf[sam]]
    if (length(wsure) == 0) 
      stop("No concordance between pred and predf on centroids data. sdifftop is too high.\n")
    coresettab <- data.frame(row.names = wsure)
    coresettab$groups <- predf[wsure]
    coresettab <- cbind(coresettab, as.matrix(tdist[wsure, 
                                                    ]))
    coresettab$disttocent <- mind[wsure]
    coresettab$difftop <- apply(tdist[wsure, ], 1, function(x) {
      diff(sort(x)[1:2])
    })
    refcoreset <- NULL
    for (g in names(tdist)) {
      L <- split(coresettab[, g], coresettab[, "groups"] == 
                   g)["TRUE"]
      inf <- lapply(L, function(x) c(median(x), max(x), 
                                     mad(x)))
      refcoreset <- cbind(refcoreset, matrix(unlist(inf), 
                                             nc = 1, dimnames = list(c("med", "max", "mad"), 
                                                                     g)))
    }
    if (is.null(sdisttocent)) {
      sdisttocent <- max(round((refcoreset[2, ] - refcoreset[1, 
                                                             ])/refcoreset[3, ]))
      
      sdisttocent <- refcoreset["med", pred] + sdisttocent * 
        refcoreset["mad", pred]
    }
    else {
      if (length(sdisttocent) == 1) 
        sdisttocent <- refcoreset["med", pred] + sdisttocent * 
          refcoreset["mad", pred]
    }
    if (verbose) {
      print(refcoreset)
      print(unique(cbind(pred, round(sdisttocent, 3))))
    }
    scoregroup <- c(`TRUE` = "CORE.OUTLIER", `FALSE` = "CORE")[as.character(mind > 
                                                                              sdisttocent)]
    scoregroup[which(difftop <= sdifftop)] <- "MIXTE"
    names(scoregroup) <- names(difftop)
    pred2 <- pred
    pred2[!scoregroup %in% "CORE"] <- NA
    if (!is.null(maxDist) & dist.meth %in% c("pearson", "spearman")) {
      pred2[which(mind > maxDist)] <- NA
      scoregroup[which(mind > maxDist)] <- "OUTLIER"
    }
    else {
      maxDist <- NULL
    }
  }
  else {
    scoregroup <- rep("ND", length(pred))
    scoregroup[which(difftop <= sdifftop)] <- "MIXTE"
    pred2 <- pred
    pred2[scoregroup %in% "MIXTE"] <- NA
    if (!is.null(maxDist) & dist.meth %in% c("pearson", "spearman")) {
      pred2[which(mind > maxDist)] <- NA
      scoregroup[which(mind > maxDist)] <- "OUTLIER"
    }
    else {
      pred2[which(mind > quantile(mind, 0.95))] <- NA
      scoregroup[which(mind > quantile(mind, 0.95))] <- "OUTLIER"
      maxDist <- quantile(mind, 0.95)
    }
  }
  tdist <- tdist[srcs, ]
  difftop <- difftop[srcs]
  mind <- mind[srcs]
  pred <- pred[srcs]
  predf <- predf[srcs]
  pred2 <- pred2[srcs]
  scoregroup <- scoregroup[srcs]
  list(dist.scores = tdist, distToNearestCentroid = mind, diffDistTopCentroids = difftop, 
       pred = pred, predwmixed = predf, predCore = pred2, group.confidence = scoregroup, 
       dist.meth = dist.meth, cutoffdiffdist = sdisttocent, 
       cutoffdisttocent = sdifftop, cutoffdistmax = maxDist)
}


cit.dist <- function (x, meth = "pearson", use = "complete.obs", diag = FALSE, 
                      upper = FALSE, p = 2, replaceNA = TRUE) 
{
  res <- NULL
  
  if (!is.na(meth) & !is.null(meth)) {
    if (meth == "propNonEqual") {
      nc <- ncol(x)
      nr <- nrow(x)
      m <- apply(x, 1, function(xi) rowSums(t(t(x) != xi), 
                                            na.rm = T) + rowSums(is.na(t(t(x) != xi))))/nc
      res <- as.dist(m, diag = diag, upper = upper)
      maxdist <- ceiling(max(res, na.rm = TRUE))
    }
    if (meth == "cosine") {
      tmp <- x %*% t(x)
      norm <- sqrt(diag(tmp)) %*% t(sqrt(diag(tmp)))
      res <- as.dist(1 - (tmp/norm), diag = diag, upper = upper)
      maxdist <- ceiling(max(res, na.rm = TRUE))
    }
    if (meth %in% c("euclidian", "maximum", "manhattan", 
                    "canberra", "binary", "minkowski")) {
      res <- dist(x, method = meth, diag = diag, upper = upper, 
                  p = p)
      maxdist <- ceiling(max(res, na.rm = TRUE))
    }
    if (meth %in% c("pearson", "spearman")) {
      maxdist <- 2
      if (use!="complete.obs") {
        use <- "all.obs"
        print("NB - function cit.dist (utile clustering) : parameter 'use' was set to default value = 'all.obs'")
      }
      res <- as.dist(1 - cor(t(x), use = use, method = meth), 
                     diag = diag, upper = upper)
    }
  }
  wNA <- which(is.na(res))
  if (length(wNA) > 0) {
    if (replaceNA) {
      res[wNA] <- maxdist
      print(paste("Warning - function cit.dist (utile clustering) : d(x,y) = NA is replaced by 0 if x=y,", 
                  maxdist, "otherwise"))
      print("NB : maybe an other distance metric (ex. euclidian, manhattan,...) could avoid inconsistencies !!")
    }
    else {
      print("Warning - function cit.dist (utile clustering) : there are NA values in the distance matrix")
    }
  }
  if (is.null(res)) 
    stop("ERROR - function cit.dist : uncorrect value for parameter 'meth' and/or 'use'")
  return(res)
}

ncc.corr<-function(centroid,newdata)
{
    ## a) reduce centroid and newdata to common genes
    genes<-intersect(rownames(centroid),rownames(newdata))
    cx<-centroid[genes,]
    bx<-newdata[genes,]
    
    ## b) calculate correlation for each centroid to each sample
    dd<-matrix(nrow=ncol(bx),ncol=ncol(cx),dimnames=list(colnames(bx),colnames(cx)))
    for (i in 1:ncol(bx)){
        for(j in 1:ncol(cx)){
            dd[i,j]<-1-cor(cx[,j],bx[,i], use = "complete.obs")
        }}
    
    #  c) find min dist, i.e. class for each sample
    pr.euc<-apply(dd,1,which.min)
    pr<-colnames(dd)[pr.euc]
    names(pr)<-names(pr.euc)
    
    return(pr)
}


options("error"=geterrmessage)
