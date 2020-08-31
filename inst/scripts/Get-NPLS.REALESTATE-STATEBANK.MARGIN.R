## Allocate the lists
Y <- list()
X <- list()

## Assign variables


Y[["real_estate"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/stock_returns/real_estate_industry_return_residuals.txt'),
                                         header = TRUE)[, "quar.st.resi.Esta", drop = FALSE])


Y[["state_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLS/state_bank_NPLs_raw_data.txt'),
                                         header = TRUE)[, "NPLs", drop = FALSE])



X[["real_estate"]] <- NA # no covariates attached for residuals
X[["state_bank"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/X/state_bank_specific_factors_data.txt'),
                                         header = TRUE)[,4:17])

X[["macroeconomics"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/X/macroeconomic_covariate_data.txt'),
                                           header = TRUE)[, 5:9])


## Post processing
lapply(Y, function(x) colnames(x) <- NULL)
