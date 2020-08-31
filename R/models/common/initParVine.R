#' Set the initial value for vine model
#'
#' Initialize the Vine parameters
#' @param Y.continuous NA
#' @param X.continuous NA
#' @param Y.discrete NA
#' @param X.discrete NA
#' @param X.FUN NA
#' @param Vine.InitTree NA
#' @param Vine.Margis.beta NA
#' @param Vine.Margis.betaIdx NA
#' @param Vine.Cpl.beta NA
#' @param Vine.Cpl.betaIdx NA
#' @param Vine.RVM NA
#' @param Mdl.MargisType NA
#' @param Mdl.betaInit NA
#' @param Mdl.parLink NA
#' @param Mdl.varSelArgs NA
#' @param Mdl.parUpdate NA
#' @param MCMC.optimInit NA
#' @return NA
#' @references Li 2017
#' @author Feng Li, Central University of Finance and Economics.
#' @export
initParVine <- function(Y.continuous = NA,
                        X.continuous = NA,
                        Y.discrete = NA,
                        X.discrete = NA,
                        X.FUN = X.FUN,
                        Vine.InitTree,
                        Vine.Margis.beta,
                        Vine.Margis.betaIdx,
                        Vine.Cpl.beta,
                        Vine.Cpl.betaIdx,
                        Vine.RVM,
                        Mdl.MargisType,
                        Mdl.betaInit,
                        Mdl.parLink,
                        Mdl.varSelArgs,
                        Mdl.parUpdate,
                        MCMC.optimInit)
{
    ## Initialize marginal parameters
    Margis.X <- list()
    Margis.Y <- list()
    Vine.Margis.parUpdate <- list()
    Vine.Cpl.parUpdate <- list()

    Margis.NM <- Vine.InitTree[["MargisNM"]]
    for(blockCaller in 1:length(Vine.Margis.beta))
    {
        blockNM <- names(Margis.NM)[blockCaller]
        Margis.X[[blockCaller]] <- list()
        Margis.Y[[blockCaller]] <- list()
        Vine.Margis.parUpdate[[blockCaller]] <- list()

        for(margiCaller in 1:length(Vine.Margis.beta[[blockCaller]]))
        {
            margiType <- Vine.InitTree[["MargisType"]][[blockCaller]][[margiCaller]]
            margiIsdiscrete <- Vine.InitTree[["isdiscrete"]][[blockCaller]][[margiCaller]]

            ## cat(blockCaller, margiCaller, "\n")
            parUpdate <- Mdl.parUpdate[margiType]
            parLink <- Mdl.parLink[margiType]
            betaInit <- Mdl.betaInit[margiType]
            varSelArgs <- Mdl.varSelArgs[margiType]

            ## cat(blockNM, blockCaller, margiCaller, "\n")

            if(margiIsdiscrete == FALSE)
            {
                X = X.continuous[[blockNM]]
                Y = Y.continuous[[blockNM]]
            }
            else
            {
                X = X.discrete[[blockNM]]
                Y = Y.discrete[[blockCaller]][ ,Margis.NM[[blockNM]][[margiCaller]],
                                              drop = FALSE]
            }
            Margis.X.curr <- rapply(X.FUN[margiType], eval.parent, n = 1,  how = "replace")
            ## Y <- Dat.Y.continuous[margiNM]

            parOut <- initPar0(Mdl.varSelArgs = Mdl.varSelArgs,
                               Mdl.betaInit = betaInit,
                               Mdl.X = Margis.X.curr,
                               Mdl.Y = NA,
                               Mdl.parLink = parLink,
                               parUpdate = parUpdate,
                               MCMC.optimInit = MCMC.optimInit)

            Margis.Y[[blockCaller]][[margiCaller]] <- Y
            Margis.X[[blockCaller]][[margiCaller]] <- Margis.X.curr[[1]] ## [[1]] to take away list name
            Vine.Margis.beta[[blockCaller]][[margiCaller]] <- parOut[["Mdl.beta"]][[margiType]]
            Vine.Margis.betaIdx[[blockCaller]][[margiCaller]] <- parOut[["Mdl.betaIdx"]][[margiType]]
            Vine.Margis.parUpdate[[blockCaller]][[margiCaller]] <- parUpdate[[1]]
        }
    }


    ## Initialize copula parameters
    family <- Vine.RVM[["family"]]
    familyVec <- CplSwapCharNum(family[lower.tri(family)])
    Cpl.X <- list()
    for(margiCaller in 1:length(Vine.Cpl.beta))
    {
        margiType <- familyVec[margiCaller]
        parUpdate <- Mdl.parUpdate[margiType]
        parLink <- Mdl.parLink[margiType]
        betaInit <- Mdl.betaInit[margiType]
        varSelArgs <- Mdl.varSelArgs[margiType]

        X <- X.continuous[[1]]
        warning("X should be from two corresponding margins")

        Cpl.X.curr <- rapply(X.FUN[margiType], eval.parent, n = 1, how = "replace")


        ## Y <- Dat.Y.continuous[margiNM]

        parOut <- initPar0(Mdl.varSelArgs = Mdl.varSelArgs,
                           Mdl.betaInit = betaInit,
                           Mdl.X = Cpl.X.curr,
                           Mdl.Y = NA,
                           Mdl.parLink = parLink,
                           parUpdate = parUpdate,
                           MCMC.optimInit = MCMC.optimInit)

        Cpl.X[[margiCaller]] <- Cpl.X.curr[[1]] ## [[1]] to take away list name
        Vine.Cpl.beta[[margiCaller]] <- parOut[["Mdl.beta"]][[margiType]]
        Vine.Cpl.betaIdx[[margiCaller]] <- parOut[["Mdl.betaIdx"]][[margiType]]
        Vine.Cpl.parUpdate[[margiCaller]] <- parUpdate[[1]]
    }

    ## browser()
    out <- list(Vine.Margis.beta = Vine.Margis.beta,
                Vine.Margis.betaIdx = Vine.Margis.betaIdx,
                Vine.Margis.parUpdate = Vine.Margis.parUpdate,
                Vine.Cpl.beta = Vine.Cpl.beta,
                Vine.Cpl.betaIdx = Vine.Cpl.betaIdx,
                Vine.Cpl.parUpdate = Vine.Cpl.parUpdate,
                Margis.Y = Margis.Y,
                Margis.X = Margis.X,
                Cpl.X = Cpl.X)
    return(out)
}


#' @export
par2RVMpar <- function(Vine.Cpl.par, family.Lst)
{
    ## Standardized from parameterization
    Vine.Cpl.parStd <-  mapply(parCplRep2Std, CplNM = family.Lst,
                               parCplRep = Vine.Cpl.par, SIMPLIFY = FALSE)
    ## Reserve space
    Vine.nPairs = length(Vine.Cpl.par)
    Vine.dim <- (1+sqrt(1+8*Vine.nPairs))/2
    nObs <- length(Vine.Cpl.par[[1]][[1]])

    par.Ary <- array(0, c(Vine.dim, Vine.dim, nObs))
    par2.Ary <- par.Ary

    mat <- matrix(0, Vine.dim, Vine.dim)
    idx <- which(lower.tri(mat))
    idx2dim <- arrayInd(idx, c(Vine.dim, Vine.dim))

    for(i in 1:Vine.nPairs)
    {
        nPar <- length(Vine.Cpl.par[[i]])
        if(nPar == 2)
        {## Only for some bivariate copula with bivariate parameters

            if(tolower(family.Lst[i]) %in% c("bb7"))
            {
                ## Special case that VineCopula parameter order is different from cdcopula
                ## for some bivariate copulas
                par.Ary[idx2dim[i, 1], idx2dim[i, 2], ] <- Vine.Cpl.parStd[[i]][[2]]
                par2.Ary[idx2dim[i, 1], idx2dim[i, 2], ] <- Vine.Cpl.parStd[[i]][[1]]
            }
            else
            {
                par.Ary[idx2dim[i, 1], idx2dim[i, 2], ] <- Vine.Cpl.parStd[[i]][[1]]
                par2.Ary[idx2dim[i, 1], idx2dim[i, 2], ] <- Vine.Cpl.parStd[[i]][[2]]
            }
        }
        else
        {
            ## One parameter copula
            par.Ary[idx2dim[i, 1], idx2dim[i, 2], ] <- Vine.Cpl.parStd[[i]][[1]]
        }
    }

    ## Assign array with std values
    out <- list(par.Ary = par.Ary, par2.Ary = par2.Ary)
    return(out)
}
