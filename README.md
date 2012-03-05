Covariate-dependent copula model
================================

Copyright 
---------

Feng Li <feng.li@stat.su.se>

Preparations
------------

### Environment requirements

*    git tools
*    R > 2.14.0
*    R compiler package (optional)

### Build a fast R with BLAS (optional)

*  See the ["R Installation and Administration"](http://cran.r-project.org/doc/manuals/R-admin.pdf) for details.

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
