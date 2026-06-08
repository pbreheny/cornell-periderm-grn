# Gene regulatory network reconstruction of zebrafish periderm

This repository contains the code and data needed to reproduce the analyses and figures from our manuscript. The workflow is managed by [Snakemake](https://snakemake.readthedocs.io), which handles all dependencies automatically: running `snakemake figs/fig3.pdf`, for example, will fit all the required network models and produce the figure.

## Quick start

To reproduce all figures:

```bash
snakemake --cores <N>
```

To reproduce a specific figure:

```bash
snakemake figs/fig3.pdf
snakemake figs/fig5a.pdf figs/fig5b.pdf figs/fig5c.csv
snakemake figs/fig6.pdf
```

Snakemake will determine which upstream steps (network fitting, etc.) are needed and run them in the correct order.

## Container (recommended)

The Snakefile is configured to run all rules inside an [Apptainer](https://apptainer.org) container, which provides a fully reproducible software environment. To build the container:

```bash
apptainer build container/periderm-grn.sif container/periderm-grn.def
```

To use the container, Snakemake must be told to enable Apptainer:

```bash
snakemake --cores <N> --software-deployment-method apptainer
```

## Without the container

If you prefer to use your own R and Python installations, you can run Snakemake without the `--software-deployment-method apptainer` flag. You will need:

**R packages** (install via `install.packages()` / `BiocManager::install()`):

| Package | Source |
|---------|--------|
| ashr | CRAN |
| data.table | CRAN |
| forcats | CRAN |
| glmnet | CRAN |
| glue | CRAN |
| ggplot2 | CRAN |
| ggraph | CRAN |
| igraph | CRAN |
| ncvreg | CRAN |
| readxl | CRAN |
| shadowtext | CRAN |
| stringr | CRAN |
| MetBrewer | CRAN |
| AnnotationDbi | Bioconductor |
| DESeq2 | Bioconductor |

**Python packages**:

```bash
pip install docopt inferelator
```

## Directory structure

```
data/
  raw/         Primary input data (expression, gene lists)
  chip/        ChIP/CUT&RUN gold-standard binding data
  go/          Gene Ontology annotations (JSON, downloaded from PANTHER)
  knockout/    Differential expression results from KO experiments
network/
  gl/          Graphical lasso networks
  inferelator/ Inferelator 3.0 networks
  wgcna/       Correlation-based (WGCNA)
figs/          Output figures (PDF) and supplementary tables (CSV)
scripts/       Analysis scripts called by Snakemake rules
R/             Shared R functions sourced by the scripts
container/     Apptainer definition file
```
