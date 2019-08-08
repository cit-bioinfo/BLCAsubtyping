## ----eval=FALSE----------------------------------------------------------
#  require(devtools)
#  devtools::install_github("cit-bioinfo/BLCAsubtyping")

## ----message=FALSE-------------------------------------------------------
library(BLCAsubtyping)
data(cit) 

## ------------------------------------------------------------------------
cl <- classify(expMat = cit$expMat, gpl = cit$gpl, symbol = "Symbol", classification.systems = c("Baylor", "UNC", "MDA", "CIT", "Lund", "TCGA"))

## ----eval=FALSE----------------------------------------------------------
#  head(cl)

## ----echo=FALSE----------------------------------------------------------
kableExtra::kable_styling(knitr::kable(head(cl), digits = 3, format = "html"))

