# BLCAsubtyping
Transcriptomic tools to classify bladder tumours according to six published molecular classifications : Baylor[1], UNC[2], MDA[3], Lund[4], CIT-Curie[5], TCGA[6]

## Citation
For now, you can cite the following bioRxiv preprint: bioRxiv 488460; doi: https://doi.org/10.1101/488460

## Consensus Class

For the Consensus BLCA classification, please use the consensusMIBC package: https://github.com/cit-bioinfo/consensusMIBC 

## Install
You may install this package with [devtools]:

[devtools]: https://github.com/hadley/devtools

```{r}
require(devtools)
devtools::install_github("cit-bioinfo/BLCAsubtyping", build_vignettes = TRUE)
```

## Documentation

After installation, please refer to the vignette or function help files by calling:

```{r}
library(BLCAsubtyping)
vignette("BLCAsubtyping")
? classify
```

## References
[1] Mo, Q. et al. Prognostic Power of a Tumor Differentiation Gene Signature for Bladder Urothelial Carcinomas. J. Natl. Cancer Inst. (2018).

[2] Damrauer, J. S. et al. Intrinsic subtypes of high-grade bladder cancer reflect the hallmarks of breast cancer biology. Proc. Natl. Acad. Sci. U.S.A. 111, 3110–3115 (2014).

[3] Choi, W. et al. Identification of distinct basal and luminal subtypes of muscle-invasive bladder cancer with different sensitivities to frontline chemotherapy. Cancer Cell 25, 152–165 (2014).

[4] Marzouka, N. et al. A validation and extended description of the Lund taxonomy for urothelial carcinoma using the TCGA cohort. Scientific Reports 8, 3737 (2018).

[5] Rebouissou, S. et al. EGFR as a potential therapeutic target for a subset of muscle-invasive bladder cancers presenting a basal-like phenotype. Sci Transl Med 6, 244ra91 (2014).

[6] Robertson, A. G. et al. Comprehensive Molecular Characterization of Muscle-Invasive Bladder Cancer. Cell 171, 540-556.e25 (2017).
 

