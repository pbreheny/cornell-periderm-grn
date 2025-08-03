# Gene regulatory network reconstruction of Zebrafish periderm

## Installation

You will need to install R, Python, and bedtools. Then

``` bash
python -m pip install inferelator
python -m pip install snakemake
python -m pip install docopt
```

and

``` bash
snakemake --cores 1 setup-r
```

which will install the necessary R packages (if they are not already installed).


## Directory organization

* `bin`: executable commands. Commands that start with `batch-` are intended to be run on Argon.
* `data`
  * `raw`: Primary data.
  * `chip`: Gold standard experiments (something of a misnomer, however, as these are actually CUT&RUN experiments)
  * `go`: Gene Ontology (GO) term analysis. json files downloaded from PANTHER.
  * `knockout`: Results from differential expression analysis of knockout experiments.
  * `network`: All the fitted networks
    * `cv`: Lasso network based on cross-validation
    * `genie`: Genie networks
    * `gl`: Graphical lasso networks (no lambda selection)
    * `inferelator`: Inferelator 3.0 networks
    * `wgcna`: Correlation networks
* `R`: R functions
* `scripts`: Intended for final-product, reproducible analyses to be shared with Rob & Sunil, published, etc. If updates break this code, they should be fixed (or demoted to `explore`, if the analysis is no longer relevant)
