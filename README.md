# Covariate-dependent copula model

## Copyright

Feng Li  <feng.li@cufe.edu.cn>
School of Statistics and Mathematics
Central University of Finance and Economics
http://feng.li/

## Preparations

### Environment requirements

* git tools (for downloading the library)

* R > 3.0.0

      * Build a fast R with BLAS (optional)

         * See the ["R Installation and Administration"](http://cran.r-project.org/doc/manuals/R-admin.pdf) for details.

      * Required R packages: `mvtnorm` `optimx` `Rmfpr` `numDeriv` `optimix`.

      * Optional R packages: `compiler`, `parallel`.

### Installations

#### Clone the repository and submodules

    git clone git@bitbucket.org:fli/cdcopula.git
    cd copula
    git submodule update --init --recursive

#### Load the package to R


* This package depends on `flutils` package ("https://bitbucket.org/fli/flutils/"). Use git to
  clone the latest version to the disk. In the path `cdcopula/R/`, make a symbolic link to
  the `flutils` package.

* Let the environment variable `R_CPL_LIB_ROOT_DIR` in `~/.Renviron` file point to the
  root directory of the the package.

### Post-installation

#### Sync with the remote repository and submodules

    cd copula
    git pull

## Run the copula model

* The example model configuration files are located at `$R_CPL_LIB_ROOT_DIR/inst/config/`.

* Running scripts are provided at `$R_CPL_LIB_ROOT_DIR/inst/bin/`.  For more information,
  you may type

        inst/bin/CplRun --help

* Within R, use

        > source("inst/bin/CplRun")

command will invoke the default setting and run the default model.

## Implement other copula models

Edit the following files

### The copula model

* `dCpl.R` Density function for copulas.

* `logCplGrad.R` Gradient function for log copula density w.r.t. copula parameters.

* `logCplRepGrad.R` Gradient function for reparameterized log copula density w.r.t. copula
parameters which may require

  * `parCplRep2Std.R`

  * `kendalltauGrad.R`

  * `lambdaGrad.R`

### The marginal model

* If you want to implement a native marginal model edit the following files

  * `MargiModel.R` CDF and PDF of the marginal distribution

  * `MargiModelGrad.R` Gradient for CDF and PDF of the marginal distribution

* If you want to include a foreign model with existing algorithms, edit the following
files


  * `MargiModelForeignEval.R`

  * `MargiModelForeignPred.R`

## References

* Li, F. 2015 _Modeling covariate-contingent correlation and tail-dependence with copulas_.

* Villani, M., Kohn, R., Nott, D., (2012) _Generalized Smooth Finite Mixtures_, Journal of Econometrics, forthcoming.
