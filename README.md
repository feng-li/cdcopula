Covariate-dependent copula model
================================

Copyright 
---------

Feng Li  <feng.li@stat.su.se>

Preparations
------------

### Environment requirements

* git tools (for downloading the library)

* R > 2.14.0
     
   * Build a fast R with BLAS (optional)

      * See the ["R Installation and Administration"](http://cran.r-project.org/doc/manuals/R-admin.pdf) for details.

   * Required R packages
      * `mvtnorm` 

   * Optional R packages
      * `compiler`
      * `parallel`

Installations
-------------

### Clone the repository and submodules
    
    git clone git@bitbucket.org:fli/copula.git
    cd copula
    git submodule update --init --recursive

### Load the package to R

    source("CplMain.R")

Post-installation
-----------------
    
* Sync with the remote repository and submodules 
    
        cd copula
        git pull
        git submodule sync 
        git submodule foreach git checkout master
        git submodule foreach git pull
        git submodule update

Run the copula model
--------------------

    CplMain(setupfile="setup.R")

    
References
----------

* Li, F. 2012 Copula paper (manuscript)
* Villani, M., Kohn, R., Nott, D., (2012) _Generalized Smooth Finite Mixtures_, Journal of Econometrics, forthcoming.
