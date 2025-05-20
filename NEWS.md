
# strpip

## v0.1.3 (20250520)

### New Features
- `df_to_list` opposite to existing `list_to_df`

### Enhancement
- Rename `import_biomart_orthologs()` as `import_biomart()`
- Rename `convert_human_to_mouse()` and `convert_human_to_mouse()` as `get_conversion_dict()`
- Rename `convert_orthologue_vector()` to `convert_genes()`
- Rename `convert_orthologue_df()` to `convert_df()`
- Rename `convert_orthologue_matrix` to `convert_exprs()`
- `write_gmt()` now takes gene list as input
- Replace `save_name` and `save_dir` argument with `file` for `write_gmt()` and `write_gct()`
- Improve efficiency of `summarize_genes()` with data.table
- Rename `inst/extdata/` files to `biodict_*` and `omnipathdb_*`

### Bugfix
- `write_gmt()` uses `writeLines()` instead of `write.table()`

## v0.1.2 (20250430)
### Bugfix
* Update `read.table(fill = TRUE)` for `read_gmt()`

## v0.1.1

### New Features
* **BiomaRt Functions**:
  - Replaced `import_biomart` with `import_biomart_mouse`, `import_biomart_human`, `import_biomart_orthologs` for wider functionality.

* **Gene Conversion Functions**:
  - Added `convert_human_to_mouse()`: Converts human gene symbols to mouse orthologs with support for various mapping methods (one-to-one, one-to-many, etc.).
  - Added `convert_mouse_to_human()`: Converts mouse gene symbols to human orthologs with similar mapping capabilities.
  - Added `convert_orthologs_vector()`: A versatile function for converting gene symbols between species, supporting one-to-one, many-to-one, and many-to-many conversions.
  - Added `convert_orthologs_matrix()`: Converts matrices of gene expression data between human and mouse, handling various mapping scenarios.
  - Made compatible with `import_biomart_orthologs`

* **Database Integration**:
  - Added release-105 database for human, mouse and orthologs for all `BiomaRt Functions`

### Documentation
* Updated pkgdown website to include new gene conversion functions and their usage examples.
* Added detailed documentation for each new function, including parameters, return values, and examples.

### Bug Fixes
* Updated `get_xy_genes` and `get_mt_genes` to retrieve all possible XY and MT genes for each organism, beyond the ortholog genes.

### Other Changes
* Improved error handling in gene conversion functions to provide clearer messages for invalid inputs.


## v0.1.1
### New features
- Initial release of the strpip package
- Added specialized gene set identification functions:
  - X/Y chromosome genes
  - Mitochondrial genes
  - Structural genes
- Gene symbol mapping between species
- Data format conversions (GCT, GMT, TSV)
- Expression data processing utilities
- Integration with OmniPathR for pathway data

### Documentation
- Added pkgdown website
- Created Quick Start Guide
- Added OmniPath Integration Guide
- Included References documentation

### Dependencies
- biomaRt for gene ID mapping
- dplyr, tidyr, tibble, stringr for data manipulation
- OmnipathR for pathway data

### Breaking changes
- None (initial release)

### Bug fixes
- None (initial release)

### Other changes
- None (initial release)

### Deprecated and defunct
- None (initial release)
