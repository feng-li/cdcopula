
## Allocate the lists
Y <- list()
X <- list()

## Assign Y variables


Y[["CH"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/China.csv'),
                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["RS"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Russia.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])
Y[["PL"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Poland.csv'),
                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["GM"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Germany.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["CZ"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Czech.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["FR"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/France.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["ES"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Spain.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])
## Y[["BL"]] = as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/MarketIntegration/Y/Belgium.csv'),
##                                        header = TRUE)[, "Returns", drop = FALSE])




## MARGINAL DISTRIBUTION COVARIATES


X[["CH"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/China.csv"),
                                            header = TRUE)[, c(1,2,4,5)])
## X[["RS"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Russia.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])
X[["PL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Poland.csv"),
                                            header = TRUE)[, c(1,2,4,5)])
## X[["GM"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Germany.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])
## X[["CZ"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Czech.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])
## X[["FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/France.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])
## X[["ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Spain.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])
## X[["BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/X/Belgium.csv"),
##                                            header = TRUE)[, c(1,2,4,5)])





###########################################################################
###
### COVARIATES FOR COPULA PARAMETERS
###
###########################################################################

## X[["covariates"]] <- NA  # NO COVARIATES ATTACHED

X[["macroeconomics"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/Common_Factors.csv"),
                                                    header = TRUE)[,c(14,15)])


## X[["CH_RS"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_RS.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_PL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_PL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_GM"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_GM.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_CZ"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_CZ.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_FR.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CH_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CH_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["RS_PL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_PL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["RS_GM"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_GM.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["RS_CZ"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_CZ.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["RS_FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_FR.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["RS_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["RS_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/RS_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["PL_GM"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/PL_GM.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["PL_CZ"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/PL_CZ.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["PL_FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/PL_FR.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["PL_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/PL_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["PL_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/PL_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["GM_CZ"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/GM_CZ.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["GM_FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/GM_FR.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["GM_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/GM_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["GM_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/GM_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["CZ_FR"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CZ_FR.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CZ_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CZ_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["CZ_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/CZ_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["FR_ES"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/FR_ES.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])
## X[["FR_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/FR_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])

## X[["ES_BL"]] <- as.matrix(read.csv(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/MarketIntegration/Covariate/ES_BL.csv"),
##                                            header = TRUE)[, c(1,2,3,5,6)])






















## Post processing
lapply(Y, function(x) colnames(x) <- NULL)
