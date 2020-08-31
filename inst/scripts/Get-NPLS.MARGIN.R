## Allocate the lists
Y <- list()
X <- list()

## Assign variables


##Y[["real_estate"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/stock_returns/real_estate_industry_return_residuals.txt'),
##                                         header = TRUE)[, "quar.st.resi.Esta", drop = FALSE])

Y[["manufracture"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/stock_returns/manufracture_industry_return_residuals.txt'),
                                               header = TRUE)[, "quar.st.resi.Manu", drop = FALSE])


##Y[["city_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/city_bank_NPLs_raw_data.txt'),
##                                     header = TRUE)[, "NPLs", drop = FALSE])

##Y[["joint_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/joint_bank_NPLs_raw_data.txt'),
##                                               header = TRUE)[, "NPLs", drop = FALSE])

Y[["state_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/state_bank_NPLs_raw_data.txt'),
                                         header = TRUE)[, "NPLs", drop = FALSE])


##X[["real_estate"]] <- NA # no covariates attached for residuals
X[["manufracture"]] <- NA


##X[["city_bank"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/NPLS_data/X/city_bank_specific_factors_data.txt"),
##                                          header = TRUE)[, 4:17])

##X[["joint_bank"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/NPLS_data/X/joint_bank_specific_factors_data.txt"),
##                                             header = TRUE)[, 4:17])

X[["state_bank"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/NPLS_data/X/state_bank_specific_factors_data.txt"),
                                         header = TRUE)[, 4:17])




X[["macroeconomics"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/NPLS_data/X/macroeconomic_covariate_data.txt"),
                                           header = TRUE)[, 5:12])


## Post processing
lapply(Y, function(x) colnames(x) <- NULL)
