# Interoperability with file formats

``` r
options(timeout = 600)
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
set.seed(123)
counts <- matrix(rnorm(100), nrow = 20, ncol = 5)
rownames(counts) <- paste0("gene", 1:20)
colnames(counts) <- paste0("sample", 1:5)
counts <- as.data.frame(counts)
counts[1:5, 1:5]
```

    ##           sample1    sample2    sample3    sample4      sample5
    ## gene1 -0.56047565 -1.0678237 -0.6947070  0.3796395  0.005764186
    ## gene2 -0.23017749 -0.2179749 -0.2079173 -0.5023235  0.385280401
    ## gene3  1.55870831 -1.0260044 -1.2653964 -0.3332074 -0.370660032
    ## gene4  0.07050839 -0.7288912  2.1689560 -1.0185754  0.644376549
    ## gene5  0.12928774 -0.6250393  1.2079620 -1.0717912 -0.220486562

``` r
dim(counts)
```

    ## [1] 20  5

## GCT files

``` r
# Save expression matrix as GCT file
temp_gct <- tempfile(fileext = ".gct")
write_gct(counts, file = temp_gct)

# Load expression matrix from GCT file
counts <- read_gct(temp_gct)
counts[1:5, 1:5]
```

    ##           sample1    sample2    sample3    sample4      sample5
    ## gene1 -0.56047565 -1.0678237 -0.6947070  0.3796395  0.005764186
    ## gene2 -0.23017749 -0.2179749 -0.2079173 -0.5023235  0.385280401
    ## gene3  1.55870831 -1.0260044 -1.2653964 -0.3332074 -0.370660032
    ## gene4  0.07050839 -0.7288912  2.1689560 -1.0185754  0.644376549
    ## gene5  0.12928774 -0.6250393  1.2079620 -1.0717912 -0.220486562

``` r
dim(counts)
```

    ## [1] 20  5

## GMT files

``` r
# Create a list of gene lists
xy.genes <- get_xy_genes(org = "mouse")
mt.genes <- get_mt_genes(org = "mouse")
vdj.genes <- get_str_genes(org = "mouse", str = c("tcr", "bcr"))
gene.list <- list(XY = xy.genes, MT = mt.genes, VDJ = vdj.genes)
```

``` r
# Save gene list as GMT file
temp_gmt <- tempfile(fileext = ".gmt")
write_gmt(gene.list, file = temp_gmt)
```

    ## GMT file written to /tmp/RtmpyXjSqY/file1cbf10bc9f6d.gmt

``` r
# Load gene list from GMT file  
gmt <- read_gmt(temp_gmt)
head(gmt)
```

    ##        XY      MT      VDJ
    ## 1 Gm28489   mt-Tf    Igll1
    ## 2 Gm21867 mt-Rnr1 Trbv13-1
    ## 3 Gm29448   mt-Tv Trbv12-2
    ## 4 Gm21127 mt-Rnr2 Trbv13-2
    ## 5 Gm28488  mt-Tl1 Trbv12-3
    ## 6 Gm20816  mt-Nd1 Trbv13-3

``` r
# Convert list to data frame
gene.df <- list_to_df(gene.list)
head(gene.df)
```

    ##        XY      MT      VDJ
    ## 1 Gm28489   mt-Tf    Igll1
    ## 2 Gm21867 mt-Rnr1 Trbv13-1
    ## 3 Gm29448   mt-Tv Trbv12-2
    ## 4 Gm21127 mt-Rnr2 Trbv13-2
    ## 5 Gm28488  mt-Tl1 Trbv12-3
    ## 6 Gm20816  mt-Nd1 Trbv13-3

``` r
# Save gene data frame as GMT file
temp_gmt <- tempfile(fileext = ".gmt")
write_gmt(gene.df, file = temp_gmt)
```

    ## GMT file written to /tmp/RtmpyXjSqY/file1cbf1c9bb569.gmt

``` r
# Load gene data frame from GMT file  
gmt <- read_gmt(temp_gmt)
head(gmt)
```

    ##        XY      MT      VDJ
    ## 1 Gm28489   mt-Tf    Igll1
    ## 2 Gm21867 mt-Rnr1 Trbv13-1
    ## 3 Gm29448   mt-Tv Trbv12-2
    ## 4 Gm21127 mt-Rnr2 Trbv13-2
    ## 5 Gm28488  mt-Tl1 Trbv12-3
    ## 6 Gm20816  mt-Nd1 Trbv13-3
