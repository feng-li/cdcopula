Covariate-dependent copula model
================================

Copyright
---------

Feng Li  <feng.li@cufe.edu.cn>

Preparations
------------

### Environment requirements

* git tools (for downloading the library)

* R > 3.0.0

      * Build a fast R with BLAS (optional)

         * See the ["R Installation and Administration"](http://cran.r-project.org/doc/manuals/R-admin.pdf) for details.

      * Required R packages: `mvtnorm` `optimx`.

      * Optional R packages: `compiler`, `parallel`.

Installations
-------------

### Clone the repository and submodules

    git clone git@bitbucket.org:fli/copula.git
    cd copula
    git submodule update --init --recursive

### Load the package to R


* This package depends on `flutils` package ("http://repo.feng.li/flutils/"). Use git to
  clone the latest version to the disk. In the `cdcopula` package, make a symbolic link to
  the `flutils` package.

* Let the environment variable `R_CPL_LIB_ROOT_DIR` in `~/.Renviron` file point to the
  root directory of the the package.

Post-installation
-----------------

### Sync with the remote repository and submodules

    cd copula
    git pull
    git submodule sync
    git submodule foreach git checkout master
    git submodule foreach git pull
    git submodule update

Run the copula model
--------------------

* The example model configuration files are located at `$R_CPL_LIB_ROOT_DIR/inst/config/`.

* Running scripts are provided at `$R_CPL_LIB_ROOT_DIR/inst/bin/`.  For more information,
  you may type

        inst/bin/CplRun --help

* Within R, use

        > source("inst/bin/CplRun")

command will invoke the default setting and run the default model.


References
----------

* Li, F. 2015 _Modeling covariate-contingent correlation and tail-dependence with copulas_.

* Villani, M., Kohn, R., Nott, D., (2012) _Generalized Smooth Finite Mixtures_, Journal of Econometrics, forthcoming.
