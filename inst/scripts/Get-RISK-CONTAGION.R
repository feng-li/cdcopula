
## Allocate the lists
Y <- list()
X <- list()

## Assign Y variables

###########################################
#
#   PUTIAN SPECIFY
#
###########################################

## Y[["002017"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_putian_information_family_company/Y/002017_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600130"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_putian_information_family_company/Y/600130_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600680"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_putian_information_family_company/Y/600680_DD_data.txt'),
##                                           header = TRUE)[, "DD", drop = FALSE])
## Y[["600776"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_putian_information_family_company/Y/600776_DD_data.txt'),
##                                           header = TRUE)[, "DD", drop = FALSE])

########################################
#
#   CET SPECIFY
#
########################################

## Y[["000014"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_electronic_techology_family_company/Y/000014_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600330"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_electronic_techology_family_company/Y/600330_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600850"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_electronic_techology_family_company/Y/600850_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600990"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'data/Credit_risk_contagion_data/China_electronic_techology_family_company/Y/600990_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])

########################################################
#
#         CEC SPECIFY
#######################################################

## Y[["000021"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/000021_DD_data.txt'),
##                                        header = TRUE)[, "DD", drop = FALSE])
## Y[["000032"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/000032_DD_data.txt'),
##                                        header = TRUE)[, "DD", drop = FALSE])
## Y[["000066"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/000066_DD_data.txt'),
##                                        header = TRUE)[, "DD", drop = FALSE])
## Y[["000727"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/000727_DD_data.txt'),
##                                         header = TRUE)[, "DD", drop = FALSE])
Y[["000748"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/000748_DD_data.txt'),
                                         header = TRUE)[, "DD", drop = FALSE])
## Y[["600171"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/600171_DD_data.txt'),
##                                         header = TRUE)[, "DD", drop = FALSE])
Y[["600536"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/600536_DD_data.txt'),
                                         header = TRUE)[, "DD", drop = FALSE])
## Y[["600755"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/600755_DD_data.txt'),
##                                          header = TRUE)[, "DD", drop = FALSE])
## Y[["600764"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/600764_DD_data.txt'),
##                                         header = TRUE)[, "DD", drop = FALSE])
## Y[["600775"]] = as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, 'inst/extdata/CreditRisk/CEI/Y/600775_DD_data.txt'),
##                                         header = TRUE)[, "DD", drop = FALSE])

## MARGINAL DISTRIBUTION COVARIATES

## PCA_SPECIFIC COVARIATES

######################################
#
#   PUTIAN
#
#####################################

## X[["002017"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/X/002017_PCA_specific_factor_data.txt"),
##                                              header = TRUE)[, 3:6])
## X[["600130"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/X/600130_PCA_specific_factor_data.txt"),
##                                              header = TRUE)[, 3:6])
## X[["600680"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/X/600680_PCA_specific_factor_data.txt"),
##                                              header = TRUE)[, 3:6])
## X[["600776"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/X/600776_PCA_specific_factor_data.txt"),
##                                              header = TRUE)[, 3:6])

##########################################
#
#         CET
#
##########################################

## X[["000014"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_electronic_techology_family_company/X/000014_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
## X[["600330"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_electronic_techology_family_company/X/600330_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
## X[["600850"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_electronic_techology_family_company/X/600850_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
## X[["600990"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_electronic_techology_family_company/X/600990_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])

###############################################
#
#      CEC
#
###############################################

## X[["000021"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/000021_PCA_specific_factor_data.txt"),
##                                            header = TRUE)[, 3:6])
## X[["000032"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/000032_PCA_specific_factor_data.txt"),
##                                            header = TRUE)[, 3:6])
## X[["000066"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/000066_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
## X[["000727"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/000727_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
X[["000748"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/000748_PCA_specific_factor_data.txt"),
                                             header = TRUE)[, 3:6])
## X[["600171"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/600171_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
X[["600536"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/600536_PCA_specific_factor_data.txt"),
                                            header = TRUE)[, 3:6])
## X[["600755"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/600755_PCA_specific_factor_data.txt"),
##                                            header = TRUE)[, 3:6])
## X[["600764"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/600764_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])
## X[["600775"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/X/600775_PCA_specific_factor_data.txt"),
##                                             header = TRUE)[, 3:6])




###########################################################################
###
### COVARIATES USED FOR COPULA PARAMETERS
###
###########################################################################

## X[["covariates"]] <- NA  #NO COVARIATES ATTACHED

## X[["macroeconomics"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/macro_common_factor_data.txt"),
##                                                    header = TRUE)[, 3:6])

#### Macroeconomic Covariates for PUTIAN
## X[["macroeconomics"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/macro_common_factor_data1.txt"),
##                                                   header = TRUE)[, 3:7])

#################################################
#
#  CEC  FOR SUPPLEMENT SAMPLE
#
#################################################

## X[["000021_000032_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_000032_covariates.txt"),
##                                                      header = TRUE)[, 3:10])
## X[["000021_000727_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_000727_covariates.txt"),
##                                                      header = TRUE)[, 3:10])
## X[["000021_000748_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_000748_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000021_600171_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_600171_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000021_600536_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_600536_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000021_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000021_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["000032_000066_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_000066_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_000727_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_000727_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_000748_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_000748_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_600171_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_600171_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_600536_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_600536_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_600755_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000032_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000032_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["000066_000727_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_000727_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000066_000748_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_000748_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000066_600171_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_600171_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000066_600536_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_600536_covariates.txt"),
##                                                      header = TRUE)[, 3:10])
## X[["000066_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000066_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["000727_000748_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000727_000748_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000727_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000727_600755_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000727_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000727_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["000748_600171_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000748_600171_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
X[["000748_600536_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000748_600536_covariates.txt"),
                                                      header = TRUE)[, 3:14])
## X[["000748_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000748_600755_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000748_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000748_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["000748_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000748_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["600171_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600171_600755_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["600171_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600171_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["600536_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600536_600755_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["600536_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600536_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])

## X[["600755_600764_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600755_600764_covariates.txt"),
##                                                      header = TRUE)[, 3:14])
## X[["600755_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600755_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:10])

## X[["600764_600775_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/600764_600775_covariates.txt"),
##                                                      header = TRUE)[, 3:10])

################################################
#
#    PUTIAN
#
###############################################

## X[["002017_600130_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/002017_600130_covariates.txt"),
##                                                     header = TRUE)[, 3:15])
## X[["002017_600680_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/002017_600680_covariates.txt"),
##                                                     header = TRUE)[, 3:15])
## X[["002017_600776_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/002017_600776_covariates.txt"),
##                                                     header = TRUE)[, 3:15])
## X[["600130_600680_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR,"data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/600130_600680_covariates.txt"),
##                                                     header = TRUE)[, 3:15])
## X[["600130_600776_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/600130_600776_covariates.txt"),
##                                                     header = TRUE)[, 3:15])
## X[["600680_600776_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "data/Credit_risk_contagion_data/China_putian_information_family_company/COVARIATES/600680_600776_covariates.txt"),
##                                                     header = TRUE)[, 3:15])

##################################################
#
#   CEC FOR MAIN STUDY SAMPLES
#
#################################################

## X[["000021_000066_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_000066_covariates.txt"),
##                                                     header = TRUE)[, 3:10])
## X[["000021_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000021_600755_covariates.txt"),
##                                                     header = TRUE)[, 3:14])
## X[["000066_600755_specific"]] <- as.matrix(read.table(file.path(CDCOPULA_LIB_ROOT_DIR, "inst/extdata/CreditRisk/CEI/COVARIATES/000066_600755_covariates.txt"),
##                                                     header = TRUE)[, 3:14])







## Post processing
lapply(Y, function(x) colnames(x) <- NULL)
