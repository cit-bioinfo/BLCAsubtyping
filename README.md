# BLCAsubtyping
Transcriptomic tools to classify bladder tumours according to six published molecular classifications : Baylor[1], UNC[2], MDA[3], Lund[4], CIT-Curie[5], TCGA[6]

## Citation
For now, please provide a link to this github repository:
<https://github.com/cit-bioinfo/BLCAsubtyping>

## Install
You may install this package with [devtools]:

[devtools]: https://github.com/hadley/devtools

```{r}
library(devtools)
devtools::install_github("cit-bioinfo/BLCAsubtyping")
library(BLCAsubtyping)
```
## Example
The package includes an example dataset [5] to illustrate the use of the main function.
```{r}
data(example_dat) 
``` 
is a list 'cit' with two items 'expMat' and 'gpl'

The function 'classify' is used to subtype a batch of transcriptomic profiles according to one or several of the 6 classification implemented.

In the following call to 'classify', the samples will be classified according to all 6 classification systems.
```{r}
cl <- classify(expMat = cit$expMat, gpl = cit$gpl, symbol = "Symbol", classification.systems = c("Baylor", "UNC", "MDA", "CIT", "Lund", "TCGA"))
#predicting Baylor subtypes......DONE 
#predicting UNC subtypes...[1] "47 of 47 genes from the initial predictor are measured in this dataset"
#123456789101112131415161718192021222324252627282930...DONE 
#predicting CIT subtypes......DONE 
#predicting Lund subtypes......DONE 
#predicting MDA subtypes......DONE 
#predicting TCGA subtypes......DONE 
``` 

The function returns a dataframe with subtyping results from each classification system for all samples.
```{r}
head(cl)
#       ID Baylor.subtype UNC.subtype CIT.subtype Lund.subtype MDA.subtype      TCGA.subtype
#1 CIT.038          Basal       Basal         MC7    Ba/Sq-Inf       basal    Basal_squamous
#2 CIT.073 Differentiated     Luminal         MC5       GU-Inf     luminal           Luminal
#3 CIT.075          Basal       Basal         MC7        Ba/Sq       basal    Basal_squamous
#4 CIT.078          Basal       Basal         MC7        Ba/Sq       basal    Basal_squamous
#5 CIT.085          Basal     Luminal         MC1    UroA-Prog     luminal Luminal_papillary
#6 CIT.100 Differentiated     Luminal         MC2       GU-Inf    p53-like           Luminal
``` 

## References
[1] Mo, Q. et al. Prognostic Power of a Tumor Differentiation Gene Signature for Bladder Urothelial Carcinomas. J. Natl. Cancer Inst. (2018).

[2] Damrauer, J. S. et al. Intrinsic subtypes of high-grade bladder cancer reflect the hallmarks of breast cancer biology. Proc. Natl. Acad. Sci. U.S.A. 111, 3110–3115 (2014).

[3] Choi, W. et al. Identification of distinct basal and luminal subtypes of muscle-invasive bladder cancer with different sensitivities to frontline chemotherapy. Cancer Cell 25, 152–165 (2014).

[4] Marzouka, N. et al. A validation and extended description of the Lund taxonomy for urothelial carcinoma using the TCGA cohort. Scientific Reports 8, 3737 (2018).

[5] Rebouissou, S. et al. EGFR as a potential therapeutic target for a subset of muscle-invasive bladder cancers presenting a basal-like phenotype. Sci Transl Med 6, 244ra91 (2014).

[6] Robertson, A. G. et al. Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell 171, 540-556.e25 (2017).
 

