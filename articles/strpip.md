# Get Started

``` r
library(strpip)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:rlang':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

    ## 
    ## Attaching package: 'data.table'

    ## The following object is masked from 'package:rlang':
    ## 
    ##     :=

    ## The following objects are masked from 'package:dplyr':
    ## 
    ##     between, first, last

    ## strpip 0.1.4

``` r
# load example expression matrix
data_file <- system.file("data", "GSE189410_expr.txt", package = "strpip")
counts <- read.table(data_file, sep = "\t", header = TRUE)
counts[1:5, 1:5]
```

    ##               BM_IgA_1 BM_IgA_2 BM_IgA_3 BM_IgA_4 BM_IgG_1
    ## 0610005C13Rik       18        8       18       11       12
    ## 0610006L08Rik        0        0       11        0        0
    ## 0610009B22Rik      102       93      124      128      126
    ## 0610009E02Rik       50       11       22       13       21
    ## 0610009L18Rik        5       10        3        4        3

``` r
dim(counts)
```

    ## [1] 53571    11

## Remove XY, MT, VDJ genes

``` r
# Remove xy-related genes from matrix
xy.genes <- get_xy_genes(org = "mouse")
head(xy.genes)
```

    ## [1] "Gm28489" "Gm21867" "Gm29448" "Gm21127" "Gm28488" "Gm20816"

``` r
counts <- counts[!rownames(counts) %in% xy.genes,]
dim(counts)
```

    ## [1] 50450    11

``` r
# Remove mt-related genes from matrix
mt.genes <- get_mt_genes(org = "mouse")
head(mt.genes)
```

    ## [1] "mt-Tf"   "mt-Rnr1" "mt-Tv"   "mt-Rnr2" "mt-Tl1"  "mt-Nd1"

``` r
counts <- counts[!rownames(counts) %in% mt.genes,]
dim(counts)
```

    ## [1] 50426    11

``` r
# Remove vdj-related genes from matrix
vdj.genes <- get_str_genes(org = "mouse", str = c("tcr", "bcr"))
head(vdj.genes)
```

    ## [1] "Igll1"    "Trbv13-1" "Trbv12-2" "Trbv13-2" "Trbv12-3" "Trbv13-3"

``` r
counts <- counts[!rownames(counts) %in% vdj.genes,]
dim(counts)
```

    ## [1] 49760    11

## Converting a vector of gene symbols

``` r
# Mouse to human
human_genes <- convert_genes(rownames(counts), org.from = "mouse", org.to = "human", one.to.many = FALSE)
head(human_genes)
```

    ## [1] "RNASEK-C17orf49" "NCBP2AS2"        "C2orf68"         "C4orf19"        
    ## [5] "C4orf54"         "C11orf58"

``` r
length(human_genes)
```

    ## [1] 18746

``` r
# Human to mouse
mouse_genes <- convert_genes(human_genes, org.from = "human", org.to = "mouse", one.to.many = FALSE)
head(mouse_genes)
```

    ## [1] "A1bg"    "A1cf"    "A2m"     "A3galt2" "A4galt"  "A4gnt"

``` r
length(mouse_genes)
```

    ## [1] 17115

## Converting a gene expression matrix

``` r
# Convert mouse gene symbols to human
human_counts <- convert_exprs(as.matrix(counts), org.from = "mouse", org.to = "human", many.to.one = TRUE)
```

    ## Removed 8 rows with NA values or invalid gene symbols

``` r
head(human_counts)[1:5, 1:5]
```

    ##         BM_IgA_3 BM_IgG_3 BM_IgG_1 BM_IgM_3 BM_IgG_2
    ## A1BG           1        0        0        0        0
    ## A1CF           0        0        0        0        0
    ## A2M            1        0        1       14        2
    ## A3GALT2        0        0        7        3        6
    ## A4GALT       116      197      103       72      242

``` r
dim(human_counts)
```

    ## [1] 17115    11

``` r
human_counts <- convert_exprs(as.matrix(counts), org.from = "mouse", org.to = "human", many.to.one = FALSE)
```

    ## Removed 8 rows with NA values or invalid gene symbols

``` r
head(human_counts)[1:5, 1:5]
```

    ##          BM_IgA_1 BM_IgA_2 BM_IgA_3 BM_IgA_4 BM_IgG_1
    ## NCBP2AS2      191      190      197      160      147
    ## C2orf68       213      275      339      284      202
    ## C4orf19         0        0        0        0        1
    ## C4orf54        10        4        9        8        2
    ## C11orf58      827      693      825      861      951

``` r
dim(human_counts)
```

    ## [1] 16116    11

## Annotate gene function from OmniPath

``` r
# run annotation requires a gene symbol column in the expression matrix
counts$gene <- rownames(counts)

# run annotation
counts.annotated <- run_annotation(counts, gene_column = "gene", org.from = "mouse", org.to = "human", one.to.many = FALSE)
```

    ## Removed 8 rows with NA values or invalid gene symbols

``` r
head(counts.annotated)
```

    ##   BM_IgA_1 BM_IgA_2 BM_IgA_3 BM_IgA_4 BM_IgG_1 BM_IgG_2 BM_IgG_3 BM_IgG_4
    ## 1       18        8       18       11       12       15       12       12
    ## 2        0        0       11        0        0        0        1        1
    ## 3      102       93      124      128      126      125      131      171
    ## 4       50       11       22       13       21       37        9       11
    ## 5        5       10        3        4        3        1        5        2
    ## 6       76       86       88      119       64       29       55       86
    ##   BM_IgM_1 BM_IgM_3 BM_IgM_4          gene mouse_chromosome human_gene_symbol
    ## 1       19       30       20 0610005C13Rik             <NA>              <NA>
    ## 2        0        1        0 0610006L08Rik             <NA>              <NA>
    ## 3      258      207      198 0610009B22Rik             <NA>              <NA>
    ## 4       22       23       19 0610009E02Rik             <NA>              <NA>
    ## 5       12       12        4 0610009L18Rik             <NA>              <NA>
    ## 6      106      118      157 0610010K14Rik               11   RNASEK-C17orf49
    ##   human_chromosome human_gene_symbol_unique human_chromosome_unique entity_type
    ## 1             <NA>                     <NA>                    <NA>        <NA>
    ## 2             <NA>                     <NA>                    <NA>        <NA>
    ## 3             <NA>                     <NA>                    <NA>        <NA>
    ## 4             <NA>                     <NA>                    <NA>        <NA>
    ## 5             <NA>                     <NA>                    <NA>        <NA>
    ## 6               17          RNASEK-C17orf49                      17     protein
    ##        category        parent consensus_score is_cs      database_cs is_tf
    ## 1          <NA>          <NA>              NA    NA             <NA>    NA
    ## 2          <NA>          <NA>              NA    NA             <NA>    NA
    ## 3          <NA>          <NA>              NA    NA             <NA>    NA
    ## 4          <NA>          <NA>              NA    NA             <NA>    NA
    ## 5          <NA>          <NA>              NA    NA             <NA>    NA
    ## 6 intracellular intracellular               1 FALSE OmniPath, ComPPI FALSE
    ##   database_tf
    ## 1        <NA>
    ## 2        <NA>
    ## 3        <NA>
    ## 4        <NA>
    ## 5        <NA>
    ## 6        <NA>

``` r
counts.annotated %>%
  .$human_gene_symbol %>%
  unique() %>%
  length()
```

    ## [1] 16032

## Session information

``` r
sessionInfo()
```

    ## R version 4.5.2 (2025-10-31)
    ## Platform: x86_64-pc-linux-gnu
    ## Running under: Ubuntu 24.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    ## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
    ##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
    ##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
    ## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
    ## 
    ## time zone: UTC
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] strpip_0.1.4        biomaRt_2.66.1      data.table_1.18.2.1
    ## [4] magrittr_2.0.4      rlang_1.1.7         stringr_1.6.0      
    ## [7] tibble_3.3.1        tidyr_1.3.2         dplyr_1.2.0        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] rappdirs_0.3.4       sass_0.4.10          generics_0.1.4      
    ##  [4] RSQLite_2.4.6        stringi_1.8.7        hms_1.1.4           
    ##  [7] digest_0.6.39        evaluate_1.0.5       fastmap_1.2.0       
    ## [10] blob_1.3.0           jsonlite_2.0.0       progress_1.2.3      
    ## [13] AnnotationDbi_1.72.0 DBI_1.3.0            httr_1.4.8          
    ## [16] purrr_1.2.1          Biostrings_2.78.0    httr2_1.2.2         
    ## [19] textshaping_1.0.5    jquerylib_0.1.4      cli_3.6.5           
    ## [22] crayon_1.5.3         XVector_0.50.0       dbplyr_2.5.2        
    ## [25] Biobase_2.70.0       bit64_4.6.0-1        withr_3.0.2         
    ## [28] cachem_1.1.0         yaml_2.3.12          tools_4.5.2         
    ## [31] memoise_2.0.1        filelock_1.0.3       BiocGenerics_0.56.0 
    ## [34] curl_7.0.0           png_0.1-8            vctrs_0.7.1         
    ## [37] R6_2.6.1             stats4_4.5.2         BiocFileCache_3.0.0 
    ## [40] lifecycle_1.0.5      Seqinfo_1.0.0        KEGGREST_1.50.0     
    ## [43] IRanges_2.44.0       S4Vectors_0.48.0     fs_1.6.7            
    ## [46] bit_4.6.0            ragg_1.5.1           pkgconfig_2.0.3     
    ## [49] desc_1.4.3           pkgdown_2.2.0        pillar_1.11.1       
    ## [52] bslib_0.10.0         glue_1.8.0           systemfonts_1.3.2   
    ## [55] xfun_0.56            tidyselect_1.2.1     knitr_1.51          
    ## [58] htmltools_0.5.9      rmarkdown_2.30       compiler_4.5.2      
    ## [61] prettyunits_1.2.0
