Covariate-dependent copula model
================================

Copyright 
---------
Feng Li <feng.li@stat.su.se>

Preparations
------------

Build a fast R with BLAS (options)

  See the ["R Installation and Administration"](http://cran.r-project.org/doc/manuals/R-admin.pdf) for details.

Installations
-------------

### Clone the repository
    
    git clone git@bitbucket.org:fli/copula.git
    cd copula
    git submodule update --init --recursive

### Update the submodules

    git submodule foreach git pull

### Load the package to R

    source("CplMain.R")

Run the copula model
--------------------

    CplMain(setupfile="setup.R")
