# Genomic representations in R

## `GRanges` class

## `GInteractions` class

## `HiCExperiment` class


```r
library(HiCExperiment)
contacts <- contacts_yeast(full = TRUE)
```

```
## snapshotDate(): 2023-02-13
```

```
## see ?HiContactsData and browseVignettes('HiContactsData') for documentation
```

```
## downloading 1 resources
```

```
## retrieving 1 resource
```

```
## loading from cache
```

```r
contacts
```

```
## `HiCExperiment` object with 8,757,906 contacts over 763 regions 
## -------
## fileName: "/root/.cache/R/ExperimentHub/4f412dfa66_7752" 
## focus: "whole genome" 
## resolutions(5): 1000 2000 4000 8000 16000
## active resolution: 16000 
## interactions: 267709 
## scores(2): count balanced 
## topologicalFeatures: compartments(0) borders(0) loops(0) viewpoints(0) centromeres(16) 
## pairsFile: N/A 
## metadata(0):
```

## Session info 


```r
sessioninfo::session_info()
```

```
## ─ Session info ───────────────────────────────────────────────────────────────
##  setting  value
##  version  R Under development (unstable) (2023-02-09 r83797)
##  os       Ubuntu 22.04.1 LTS
##  system   x86_64, linux-gnu
##  ui       X11
##  language (EN)
##  collate  en_US.UTF-8
##  ctype    en_US.UTF-8
##  tz       America/New_York
##  date     2023-02-13
##  pandoc   2.19.2 @ /usr/local/bin/ (via rmarkdown)
## 
## ─ Packages ───────────────────────────────────────────────────────────────────
##  package                * version   date (UTC) lib source
##  AnnotationDbi            1.61.0    2022-11-01 [1] Bioconductor
##  AnnotationHub          * 3.7.1     2023-02-06 [1] Bioconductor
##  assertthat               0.2.1     2019-03-21 [1] CRAN (R 4.3.0)
##  Biobase                  2.59.0    2022-11-01 [1] Bioconductor
##  BiocFileCache          * 2.7.1     2022-12-09 [1] Bioconductor
##  BiocGenerics           * 0.45.0    2022-11-01 [1] Bioconductor
##  BiocIO                   1.9.2     2023-01-19 [1] Bioconductor
##  BiocManager              1.30.19   2022-10-25 [1] CRAN (R 4.3.0)
##  BiocParallel             1.33.9    2022-12-23 [1] Bioconductor
##  BiocVersion              3.17.1    2022-11-04 [1] Bioconductor
##  Biostrings               2.67.0    2022-11-01 [1] Bioconductor
##  bit                      4.0.5     2022-11-15 [1] CRAN (R 4.3.0)
##  bit64                    4.0.5     2020-08-30 [1] CRAN (R 4.3.0)
##  bitops                   1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
##  blob                     1.2.3     2022-04-10 [1] CRAN (R 4.3.0)
##  bookdown                 0.32      2023-01-17 [1] CRAN (R 4.3.0)
##  bslib                    0.4.2     2022-12-16 [1] CRAN (R 4.3.0)
##  cachem                   1.0.6     2021-08-19 [1] CRAN (R 4.3.0)
##  cli                      3.6.0     2023-01-09 [1] CRAN (R 4.3.0)
##  codetools                0.2-19    2023-02-01 [2] CRAN (R 4.3.0)
##  crayon                   1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
##  curl                     5.0.0     2023-01-12 [1] CRAN (R 4.3.0)
##  DBI                      1.1.3     2022-06-18 [1] CRAN (R 4.3.0)
##  dbplyr                 * 2.3.0     2023-01-16 [1] CRAN (R 4.3.0)
##  DelayedArray             0.25.0    2022-11-01 [1] Bioconductor
##  digest                   0.6.31    2022-12-11 [1] CRAN (R 4.3.0)
##  dplyr                    1.1.0     2023-01-29 [1] CRAN (R 4.3.0)
##  ellipsis                 0.3.2     2021-04-29 [1] CRAN (R 4.3.0)
##  evaluate                 0.20      2023-01-17 [1] CRAN (R 4.3.0)
##  ExperimentHub          * 2.7.0     2022-11-01 [1] Bioconductor
##  fansi                    1.0.4     2023-01-22 [1] CRAN (R 4.3.0)
##  fastmap                  1.1.0     2021-01-25 [1] CRAN (R 4.3.0)
##  filelock                 1.0.2     2018-10-05 [1] CRAN (R 4.3.0)
##  generics                 0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
##  GenomeInfoDb             1.35.15   2023-02-02 [1] Bioconductor
##  GenomeInfoDbData         1.2.9     2023-02-13 [1] Bioconductor
##  GenomicRanges            1.51.4    2022-12-15 [1] Bioconductor
##  glue                     1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
##  HiCExperiment          * 0.99.9    2023-02-13 [1] Github (js2264/HiCExperiment@7f73f79)
##  HiContactsData         * 1.1.9     2023-02-13 [1] Github (js2264/HiContactsData@275fab4)
##  htmltools                0.5.4     2022-12-07 [1] CRAN (R 4.3.0)
##  httpuv                   1.6.8     2023-01-12 [1] CRAN (R 4.3.0)
##  httr                     1.4.4     2022-08-17 [1] CRAN (R 4.3.0)
##  InteractionSet           1.27.0    2022-11-01 [1] Bioconductor
##  interactiveDisplayBase   1.37.0    2022-11-01 [1] Bioconductor
##  IRanges                  2.33.0    2022-11-01 [1] Bioconductor
##  jquerylib                0.1.4     2021-04-26 [1] CRAN (R 4.3.0)
##  jsonlite                 1.8.4     2022-12-06 [1] CRAN (R 4.3.0)
##  KEGGREST                 1.39.0    2022-11-01 [1] Bioconductor
##  knitr                    1.42      2023-01-25 [1] CRAN (R 4.3.0)
##  later                    1.3.0     2021-08-18 [1] CRAN (R 4.3.0)
##  lattice                  0.20-45   2021-09-22 [2] CRAN (R 4.3.0)
##  lifecycle                1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
##  magrittr                 2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
##  Matrix                   1.5-3     2022-11-11 [2] CRAN (R 4.3.0)
##  MatrixGenerics           1.11.0    2022-11-01 [1] Bioconductor
##  matrixStats              0.63.0    2022-11-18 [1] CRAN (R 4.3.0)
##  memoise                  2.0.1     2021-11-26 [1] CRAN (R 4.3.0)
##  mime                     0.12      2021-09-28 [1] CRAN (R 4.3.0)
##  pillar                   1.8.1     2022-08-19 [1] CRAN (R 4.3.0)
##  pkgconfig                2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
##  png                      0.1-8     2022-11-29 [1] CRAN (R 4.3.0)
##  promises                 1.2.0.1   2021-02-11 [1] CRAN (R 4.3.0)
##  purrr                    1.0.1     2023-01-10 [1] CRAN (R 4.3.0)
##  R6                       2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
##  rappdirs                 0.3.3     2021-01-31 [1] CRAN (R 4.3.0)
##  Rcpp                     1.0.10    2023-01-22 [1] CRAN (R 4.3.0)
##  RCurl                    1.98-1.10 2023-01-27 [1] CRAN (R 4.3.0)
##  rhdf5                    2.43.0    2022-11-01 [1] Bioconductor
##  rhdf5filters             1.11.0    2022-11-01 [1] Bioconductor
##  Rhdf5lib                 1.21.0    2022-11-01 [1] Bioconductor
##  rlang                    1.0.6     2022-09-24 [1] CRAN (R 4.3.0)
##  rmarkdown                2.20      2023-01-19 [1] CRAN (R 4.3.0)
##  RSQLite                  2.2.20    2022-12-22 [1] CRAN (R 4.3.0)
##  S4Vectors                0.37.3    2022-12-07 [1] Bioconductor
##  sass                     0.4.5     2023-01-24 [1] CRAN (R 4.3.0)
##  sessioninfo              1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
##  shiny                    1.7.4     2022-12-15 [1] CRAN (R 4.3.0)
##  strawr                   0.0.9     2021-09-13 [1] CRAN (R 4.3.0)
##  SummarizedExperiment     1.29.1    2022-11-04 [1] Bioconductor
##  tibble                   3.1.8     2022-07-22 [1] CRAN (R 4.3.0)
##  tidyselect               1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
##  tzdb                     0.3.0     2022-03-28 [1] CRAN (R 4.3.0)
##  utf8                     1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
##  vctrs                    0.5.2     2023-01-23 [1] CRAN (R 4.3.0)
##  vroom                    1.6.1     2023-01-22 [1] CRAN (R 4.3.0)
##  withr                    2.5.0     2022-03-03 [1] CRAN (R 4.3.0)
##  xfun                     0.37      2023-01-31 [1] CRAN (R 4.3.0)
##  xtable                   1.8-4     2019-04-21 [1] CRAN (R 4.3.0)
##  XVector                  0.39.0    2022-11-01 [1] Bioconductor
##  yaml                     2.3.7     2023-01-23 [1] CRAN (R 4.3.0)
##  zlibbioc                 1.45.0    2022-11-01 [1] Bioconductor
## 
##  [1] /usr/local/lib/R/site-library
##  [2] /usr/local/lib/R/library
## 
## ──────────────────────────────────────────────────────────────────────────────
```
