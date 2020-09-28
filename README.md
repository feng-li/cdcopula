# Covariate-dependent copula models

## Copyright

    Feng Li  <feng.li@cufe.edu.cn>
    School of Statistics and Mathematics
    Central University of Finance and Economics
    http://feng.li/

## Cite the package

* Li, F., and Kang, Y.(2018). [Improving forecasting performance using covariate-dependent
  copula models](http://dx.doi.org/10.1016/j.ijforecast.2018.01.007), International
  Journal of Forecasting, 34(3), pp. 456-476.

## Prerequisites

### Environment requirements

* git tools (for downloading the library)

* R > 3.0.0 with the following required packages:

    * `devtools`

    * `mvtnorm` `Rmfpr` `numDeriv` `optimix`

    * `rugarch` `rmgarch`  `fGarch` `stochvol` `teigen`

    * `VineCopula`

    * `snow`, `Rmpi` (optional)

### Installation

The package is currently under developments. We are not yet ready to make it as a standard
R package. Please follow the following instructions to "load" the functionality.

#### Clone the repository and its dependencies

- This package depends on [`flutils`](https://github.com/feng-li/flutils) package. Install it with `devtools`

``` r
devtools::install_github("feng-li/flutils")
```

``` r
devtools::install_github("feng-li/cdcopula")
```


### Run a covariate-dependent copula model

* The example model configuration files are located at `inst/config/`.

* Linux executable scripts are provided at `inst/bin/`.  For more information, you may
  issue the following command in a terminal

``` sh
inst/bin/CplRun --help
```

* Or within R, use

``` R
library("cdcopula")
CplRun(Mdl.ConfigFile=file.path(system.file(package = "cdcopula"), 'config/config.BB7.GARCH.SP100-SP600.R'))
```

### Cluster Parallelization

* Run with SLURM cluster

``` sh
mpirun -np 1  inst/bin/CplRun 4 inst/config/${CONFIG_FILE}
```

And submit to slurm with

``` sh
sbatch inst/config/slurm.sh
```

## Implementing a new covariate-dependent copula model

First you may need to consult the reference paper. Then you could edit the following files
accordingly.

### for the copula model

* `pCpl.R`, `dCpl.R` Copula functions and their densities.

* `logCplGrad.R` Gradient function for log copula density w.r.t. copula parameters.

* `logCplRepGrad.R` Gradient function for reparameterized log copula density w.r.t. copula
parameters which may require

    * `parCplRep2Std.R`

    * `kendalltau.R`, `kendalltauGrad.R`

    * `lambda.R`, `lambdaGrad.R`

### for the marginal model

* If you want to implement a native marginal model edit the following files

    * `MargiModel.R` CDF and PDF of the marginal distribution

    * `MargiModelGrad.R` Gradient for CDF and PDF of the marginal distribution

    * `MargiModelInv.R` The inverse for marginal model

* If you want to include a foreign marginal model with existing algorithms, edit the
following files

    * `MargiModelForeignEval.R`

    * `MargiModelForeignPred.R`

    * `logCplPredict.R`

* If you want to include a foreign multivariate model for model comparison, edit the
  following files

    * `ModelForeignEval.R`

    * `ModelForeignPred.R`

## References

* Li, F., Panagiotelis, A. and Kang, Y. (2019) _Modelling multivariate tail-dependence
  with covariate-dependent vine copulas_, Working paper.

* Li, F., and Kang, Y.(2018). [Improving forecasting performance using covariate-dependent
  copula models](http://dx.doi.org/10.1016/j.ijforecast.2018.01.007), International
  Journal of Forecasting, 34(3), pp. 456-476.


## Acknowledgments

Feng Li's research were supported by the National Natural Science Foundation of China
(No. 11501587).
