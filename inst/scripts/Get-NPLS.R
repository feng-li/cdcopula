## Allocate the lists
Y <- list()
X <- list()

## Assign variables


##Y[["real_estate"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/stock_returns/real_estate_industry_return_residuals.txt'),
##                                         header = TRUE)[, "quar.st.resi.Esta", drop = FALSE])


Y[["manufracture"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/stock_returns/manufracture_industry_return_residuals.txt'),
                                              header = TRUE)[, "quar.st.resi.Manu", drop = FALSE])

##Y[["construct"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/construct_industry_return_residuals.txt'),
##                                         header = TRUE)[, "quar.st.resi.Const", drop = FALSE])
##Y[["finance"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/financial_industry_return_residuals.txt'),
##                                           header = TRUE)[, "quar.st.resi.Finan", drop = FALSE])
##Y[["information_technology"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/information_technology_industry_return_residuals.txt'),
##                                                         header = TRUE)[, "quar.st.resi.IT", drop = FALSE])
##Y[["manufracture"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/manufracturing_industry_return_residuals.txt'),
##                                                 header = TRUE)[, "quar.st.resi.Manu", drop = FALSE])
##Y[["mining"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/mining_industry_return_residuals.txt'),
##                                         header = TRUE)[, "quar.st.resi.Min", drop = FALSE])
##Y[["transport_warehouse"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS/Y/stock_returns/transport_and_warehousing_industry_return_residuals.txt'),
##                                           header = TRUE)[, "quar.st.resi.TW", drop = FALSE])
##Y[["wholesale_retail"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIT, 'data/NPLS/Y/stock_returns/wholesale_and_retailing_industry_return_residuals.txt'),
##                                       header = TRUE)[, "quar.st.resi.WR", drop = FALSE])


##Y[["city_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/city_bank_NPLs_residuals_data.txt'),
##                                            header =TRUE)[, "NPLs.regress.residuals", drop = FALSE])

Y[["state_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/state_bank_NPLs_residuals_data.txt'),
                                            header = TRUE)[, "NPLs.regress.residuals", drop = FALSE])

##Y[["joint_bank"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/NPLS_data/Y/NPLs/joint_bank_NPLs_residuals_data.txt'),
##                                             header = TRUE)[, "NPLs.regress.residuals", drop = FALSE])


##X[["real_estate"]] <- NA # no covariates attached for residuals
X[["manufracture"]] <- NA
##X[["city_bank"]] <- NA
##X[["joint_bank"]] <- NA
X[["state_bank"]] <- NA

X[["macroeconomics"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/NPLS_data/X/macroeconomic_covariate_data.txt"),header = TRUE)[, 5:12])

## Post processing
lapply(Y, function(x) colnames(x) <- NULL)
