## Boot
# devtools::document()
# roxygen2::roxygenise()
# devtools::check()
# devtools::install()
# devtools::build()
# usethis::use_git()
# pkgdown::build_site()
# usethis::use_github_action()
# usethis::use_pkgdown_github_pages()
# usethis::use_bioc_github_action()


## one time
# pak::pak("r-lib/httr2")
# https://stackoverflow.com/questions/72189273/r-cmd-check-github-actions-workflow-failing-on-warnings-notes
# - uses: r-lib/actions/check-r-package@v1
#        build_args: 'c("--no-manual", "--no-build-vignettes")'
#        with:
#          error-on: '"error"'


#---
#title: Database
#vignette: >
#  %\VignetteIndexEntry{Interoperability with file formats}
#  %\VignetteEngine{knitr::rmarkdown}
#  %\VignetteEncoding{UTF-8}
#author: "`r 'Author: Matthew Hung'`"
#date: "`r paste('Last Updated:', format(Sys.Date(), '%Y/%m/%d'))`"
#output:
#  html_document:
#    code_folding: hide
#knitr:
#  opts_chunk:
#    collapse: true
#    comment: '#>'
#---