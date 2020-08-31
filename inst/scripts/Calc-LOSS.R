#! /usr/bin/env Rscript
## This script calculates the MSE for tail-dependence with DGP DATA

args <- commandArgs(TRUE)

if(interactive())
{
    rm(list = ls())
    ModelFileArray <- "~/running/JOB4828.M1M2+SPLITTSPLITTBB7+nObs1000nCross4MCMC.nIter1000joint+20170825@23.15.a532ae.Rdata"

    lossfun <- "squared"
    ModelOutFile <- "~/running/old/"

}else
{
    ModelIDFile <- args[1]

    # ModelIDFile = "~/running/LOSS-JOBID.csv"


    lossfun <- args[2]

    ModelOutFile <- args[3]

    runningfiles <- dir(ModelOutFile)

    ModelIDArray <- as.matrix(read.table(ModelIDFile, header = TRUE))


    ModelFileArray <- array(NA, dim(ModelIDArray))
    for(i in 1:length(ModelIDArray))
    {
        if(!is.na(ModelIDArray[i]))
            {
                ModelFileArray[i] <- runningfiles[grep(
                    paste("JOB", ModelIDArray[i], ".*.Rdata", sep = ""), runningfiles)]
            }
    }


}

## In-sample LOSS or Out-of-sample LOSS
testing <- "testing"
##testing <- "training"

## lossfun <- "squared"
## lossfun <- "absolute"

LOSS.MAT <- array(NA, c(10, length(ModelIDArray)))

whichModel <- 0
for(Model in ModelFileArray)
{
    whichModel <- whichModel+1

    if(!is.na(Model))
    {
        ## Load the results

        load(file.path(ModelOutFile, Model))

        ## Load the code
        sourceDir("~/code/cdcopula/R", recursive = TRUE)


        ## load model config
        Mdl.ConfigEnv <- OUT.FITTED[[1]][["Mdl.ConfigEnv"]]

        if(length(Mdl.ConfigEnv) == 0)
        {
            ## Old model structure
            MCMC.Update <-OUT.FITTED[[1]][["MCMC.Update"]]
            MCMC.nIter <-OUT.FITTED[[1]][["MCMC.nIter"]]
            MCMC.burninProp <-OUT.FITTED[[1]][["MCMC.burninProp"]]
            Mdl.X <-OUT.FITTED[[1]][["Mdl.X"]]
            Mdl.parLink =OUT.FITTED[[1]][["Mdl.parLink"]]
            Mdl.crossValidIdx <-OUT.FITTED[[1]][["Mdl.crossValidIdx"]]
            nCross <-OUT.FITTED[[1]][["nCross"]]
            MdlDGP.par <-OUT.FITTED[[1]][["MdlDGP.par"]]


        }else
        {
            ## New
            MCMC.Update <- Mdl.ConfigEnv[["MCMC.Update"]]
            MCMC.nIter <- Mdl.ConfigEnv[["MCMC.nIter"]]
            MCMC.burninProp <- Mdl.ConfigEnv[["MCMC.burninProp"]]
            Mdl.X <- Mdl.ConfigEnv[["Mdl.X"]]
            Mdl.parLink = Mdl.ConfigEnv[["Mdl.parLink"]]
            Mdl.crossValidIdx <- Mdl.ConfigEnv[["Mdl.crossValidIdx"]]
            nCross <- Mdl.ConfigEnv[["nCross"]]
            MdlDGP.par <- Mdl.ConfigEnv[["MdlDGP.par"]]

        }

        if(length(MdlDGP.par) == 0)
        {
            stop("Predictive loss for parameters is only available for DGP data!")
        }



        LOSS <- list()
        for (iCross in 1:nCross)
        {

            MCMC.beta <- OUT.FITTED[[iCross]][["MCMC.beta"]]
            MCMC.sampleIdx = (MCMC.nIter*MCMC.burninProp+1):MCMC.nIter


            Mdl.Idx.testing <- Mdl.crossValidIdx[[testing]][[iCross]]


            subsetFun <- function(x, idx) {x[idx, , drop = FALSE]}
            Mdl.X.testing <- rapply(object=Mdl.X, f = subsetFun,
                                    idx = Mdl.Idx.testing, how = "replace")

            MCMC.par.testing <- parCplMCMC(MCMC.beta = MCMC.beta,
                                           Mdl.X = Mdl.X.testing,
                                           Mdl.parLink = Mdl.parLink,
                                           MCMC.Update = MCMC.Update,
                                           MCMC.sampleIdx = MCMC.sampleIdx) # use full chain
            ## browser()
            MdlDGP.par.testing <- rapply(MdlDGP.par, function(x) x[Mdl.Idx.testing, drop = FALSE], how = "replace")

            loss <- MCMC.Update
            for(CompCaller in names(MCMC.Update))
            {
                for(parCaller in names(MCMC.Update[[CompCaller]]))
                {

                    if(lossfun == "squared")
                        {
                            loss[[CompCaller]][[parCaller]] <-
                                (
                                    MCMC.par.testing[[CompCaller]][[parCaller]][, , 1]-
                                    matrix(MdlDGP.par.testing[[CompCaller]][[parCaller]],
                                           length(MCMC.sampleIdx), length(Mdl.Idx.testing),
                                   byrow = TRUE)

                        )^2
                        }

                    else if(lossfun == "absolute")
                    {
                        loss[[CompCaller]][[parCaller]] <-
                            abs(
                                MCMC.par.testing[[CompCaller]][[parCaller]][, , 1]-
                                matrix(MdlDGP.par.testing[[CompCaller]][[parCaller]],
                                       length(MCMC.sampleIdx), length(Mdl.Idx.testing),
                                       byrow = TRUE)
                            )

                    }
                    else
                    {
                        stop("No such loss function")
                    }

                }
            }
            LOSS[[iCross]] = rapply(loss, mean, how = "replace")
        }
    }
    else
    {
        LOSS <- NA
    }


    OUT.PRED.LOSS4PAR <- matrix(apply(matrix(unlist(LOSS), ncol = 4), 1, mean))

    if(length(ModelFileArray) == 1)
    {
        rownames(OUT.PRED.LOSS4PAR) <- names(unlist(MCMC.Update))
        colnames(OUT.PRED.LOSS4PAR) <- "MSE LOSS"

        cat("\nMean Square Error LOSS ",
            ifelse(testing == "testing", "(out-of-sample)", "(in-sample)"),
            "\nfor copula features:\n", sep = "")
        print(OUT.PRED.LOSS4PAR)
    }

    LOSS.MAT[, whichModel] <- OUT.PRED.LOSS4PAR
    cat(Model, "......", whichModel, "/", length(ModelFileArray),  " done\n", sep = "")

}

LOSS.ARY <- array(LOSS.MAT, c(10, dim(ModelFileArray)))
filestr <- paste("~/running/LOSS/LOSS-(", date(), ")", sep = "")

save.image(file = paste(filestr, ".Rdata", sep = ""))

# load("~/LOSS-ARY(Sun Sep 10 16:11:49 2017).Rdata")

pdf(file = paste(filestr, "." ,lossfun, ".pdf", sep = ""), width = 8, height = 3)
# pdf(file = "~/LOSS-absolute-mean.pdf", width = 8, height = 3)
# pdf(file = "~/LOSS.pdf", width = 8, height = 3)

par(mfrow = c(2, 5), mar = c(2, 1, 1, 1))

for(i in c(1:4, 9, 5:8, 10)) boxplot(LOSS.ARY[i, , ])

dev.off()
