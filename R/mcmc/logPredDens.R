logPredDens <- function(Y, x, logpost.fun.name, crossvaid.struc, splineArgs, priorArgs,
                 OUT.Params, Params_Transform, burn.in, LPDS.sampleProp)
  {
    n <- dim(Y)[1] # no. of obs.
    nIter <- dim(OUT.Params[[1]])[3] # no. of iterations
    n.burnin <- floor(nIter*burnin)
    nUsed <- nIter - n.burnin
    nCross <- length(crossvalid.struc[["training"]])

    ## The sample indices for LPDS after burn-in Make a decreasing sequence and
    ## sort it by increasing order. Just to make sure last draw is allays used
    ## The exact size may not same as the result from sample proportion.

    LPDS.sampleIdx <- sort(seq(nIter, (nIter-nUsed+1), by = -round(1/LPDS.sampleProp)))
    nSample <- length(LPDS.sampleIdx)


    if(nCross  == 1)
      {
        LPDS  = NA
        nseLPDS = NA
        logPredMatrix <- NA
      }
    else if(nCross  != 1)
      {

        ## Detect how many folds used
        if(length(crossvalid.struc[["training"]][[nCross]]) +
           length(crossvalid.struc[["testing"]][[nCross]])  == n)
          { nFold <- nCross}
        else
          { nFold <- nCross -1}

        ## The log predictive matrix
        logPredMatrix <- matrix(NA, nSample, nFold)

###----------------------------------------------------------------------------
### Calculate the predictive densities for all folds
###----------------------------------------------------------------------------
        for(iCross in 1:nFold)
          {
                                        # iTraining <- crossvalid.struc[["training"]][[iCross]]
            iTesting <- crossvalid.struc[["testing"]][[iCross]]

            ## Y.iTraining <- Y[iTraining, , drop = FALSE]
            ## x.iTraining <- x[iTraining, , drop = FALSE]

            Y.iTesting <- Y[iTesting, , drop = FALSE]
            x.iTesting <- x[iTesting, , drop = FALSE]

            which.j <- 0
            for(j in LPDS.sampleIdx) ## Just the likelihood function with posterior samples
              {

                Params.j <- lapply(OUT.Params, function(x) apply(x[, , j, iCross, drop =
                                                                   FALSE], c(1, 2), "["))
                caller.log.like <- call(logpost.fun.name,Y = Y.iTesting, x = x.iTesting,
                                        Params = Params.j, callParam = list(id =
                                                             "likelihood"), priorArgs =
                                        priorArgs, splineArgs = splineArgs, Params_Transform
                                        = Params_Transform)
                log.like <- eval(caller.log.like)

                logPost(CplNM, Mdl.Y, Mdl.X,Mdl.beta,Mdl.betaIdx,Mdl.parLink,
                        varSelArgs,MargisTypes,priArgs,parUpdate,
                        staticArgs, staticArgsOnly = FALSE, parUpdate4Pri = parUpdate)


                logPost.propOut <- logPost(
                    CplNM = CplNM,
                    Mdl.Y = Mdl.Y,
                    Mdl.X = Mdl.X,
                    Mdl.beta = Mdl.beta.prop,
                    Mdl.betaIdx = Mdl.betaIdx.prop,
                    Mdl.parLink = Mdl.parLink,
                    varSelArgs = varSelArgs,
                    MargisTypes = MargisTypes,
                    priArgs = priArgs,
                    parUpdate = parUpdate,
                    staticArgs = staticArgs)

                which.j <- which.j + 1
                logPredMatrix[which.j, iCross] <- log.like

                ## Simple progress bar
                ## progressbar(((iCross-1)*nSample + which.j), nFold*nSample)
              }

          }
      }
  }
