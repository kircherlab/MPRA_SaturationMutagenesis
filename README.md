# Saturation mutagenesis MPRA data access portal

Shiny webapp of [https://mpra.gs.washington.edu/satMutMPRA](https://mpra.gs.washington.edu/satMutMPRA) or [https://kircherlab.bihealth.org/satMutMPRA/](https://kircherlab.bihealth.org/satMutMPRA/).

## Installation

### Option 1: conda

If not already done install [miniconda](https://docs.conda.io/en/latest/miniconda.html)

set up the channels:

```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

Then run `conda install mpra-data-access-portal`

The latest version on conda is `0.1.9` so you can run `conda install mpra-data-access-portal=0.1.9` to install exactly this version.

now you can run the data portal via the command: `mpra-data-access-portal`

The default port is `8080` and the default host `0.0.0.0`. If you want to change this please change the two variables `SHINY_PORT` and export `SHINY_HOST`. E.g.:

```
export SHINY_PORT=7001
export SHINY_HOST=127.0.0.1
```

Then run `mpra-data-access-portal` again.

### Option 2: from source

install dependencies:

R version 4.3.3

R packages:
- dplyr version 1.1.4
- DT version 0.33
- ggplot2 version 3.5.1
- shiny version 1.9.1
- shinytest2 version 0.3.1
- shinyvalidate 0.1.3
- htmlwidgets version 1.6.4
- readr version 2.1.5
- stringr version 1.5.1
- plotly version 4.10.4

You can install also via mamba/conda using: `mamba env create -f environment.yaml -n mpra-data-access-portal`


Get the latest version via version release on github: 
https://github.com/kircherlab/MPRA_SaturationMutagenesis/archive/v0.1.9.tar.gz 
or cloning the repository: 
https://github.com/kircherlab/MPRA_SaturationMutagenesis/archive/v0.1.9.tar.gz and change to the version tag `v0.1.9`.


the goto the directory and run the shiny server (here on port `8080` and on host `0.0.0.0`, if not alredy defined:


```
export SHINY_HOST=${SHINY_HOST-0.0.0.0}
export SHINY_PORT=${SHINY_PORT-8080}

Rscript -e "shiny::runApp(port = ${SHINY_PORT}, launch.browser = FALSE, host = '${SHINY_HOST}')"
```
