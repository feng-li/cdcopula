#' Link function constrains for copula models
#'
#' Constrains of link functions
#' @param Mdl.X  NA
#' @param Mdl.parLink NA
#' @param Mdl.beta  NA
#' @param parUpdate  NA
#' @param Mdl.par  NA
#' @return NA
#' @references Li 2012
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @export
parCplMeanFun <- function(Mdl.X,  Mdl.parLink, Mdl.beta, parUpdate, Mdl.par)
{
    ## BB7 copula (tau, lambdaL) representation requires that
    ## 0 < lambdaL < 2^(1/2-1/(2*tau)) See Li 2012 for the proof.
    ## The condition should be calculated in the end for safety.

    ## Mdl.MargisNM <- names(Mdl.beta)
    ## if(length(Mdl.MargisNM) == 0)
    ## {
    ##     stop("name must be assigned for Mdl.beta")
    ## }

    ## CompCpl <- Mdl.MargisNM[length(Mdl.MargisNM)]

    ## Update all parameter component if missing parUpdate
    if(missing(parUpdate))
    {
        parUpdate <- rapply(Mdl.beta, function(x) TRUE, how = "replace")
    }
    if(missing(Mdl.par))
    {
        Mdl.par <- rapply(Mdl.beta, function(x) NA, how = "replace")
    }

    ## Which parameter are conditionally considered Hard coded, maybe should treat it as
    ## an input condPar <- "tau"
    condPar <- NULL

###----------------------------------------------------------------------------
### (1) update all the independent linkages
###----------------------------------------------------------------------------

    for(CompCaller in names(Mdl.beta))
    {
        ## parNM <- names(parUpdate[[CompCaller]])

        ## nPar <- length(parUpdate[[CompCaller]])
        ## updateIdx <- (unlist(parUpdate[[CompCaller]]) == TRUE)
        for(parCaller in names(Mdl.beta[[CompCaller]]))
        {
            ## Check if particular constrain is needed.
            ## if(!(tolower(parCaller) %in% tolower(condPar)))
            ## {
            ## Update the parameters for the updated part
            parCurr <- parMeanFun(X = Mdl.X[[CompCaller]][[parCaller]],
                                  beta = Mdl.beta[[CompCaller]][[parCaller]],
                                  linkArgs = Mdl.parLink[[CompCaller]][[parCaller]])

            Mdl.par[[CompCaller]][[parCaller]] <- parCurr
            ## }
        }
    }
    ## browser()
###----------------------------------------------------------------------------
### (2) Special case for conditional linkage
###----------------------------------------------------------------------------
    ## CplNM0 <- strsplit(CplNM, split = "_")[[1]][1]
    ## if(tolower(CplNM0) == "bb7" && length(condPar) != 0)
    ## {
    ##     ## The parameter tau is updated individually. NOTE that the parameter
    ##     ## tau depends on lambdaL. So when lambdal is updated, the information
    ##     ## of tau should also be updated.
    ##     if(parUpdate[[CompCpl]][["tau"]] == TRUE |
    ##        parUpdate[[CompCpl]][["lambdaL"]] == TRUE)
    ##     {
    ##         ## linkCurr <- Mdl.parLink[[CompCpl]][["lambdaL"]]
    ##         ## XCurr <- Mdl.X[[CompCpl]][["lambdaL"]]
    ##         ## betaCurr <- Mdl.beta[[CompCpl]][["lambdaL"]]

    ##         linkCurr <- Mdl.parLink[[CompCpl]][["tau"]]
    ##         XCurr <- Mdl.X[[CompCpl]][["tau"]]
    ##         betaCurr <- Mdl.beta[[CompCpl]][["tau"]]

    ##         if(tolower(linkCurr[["type"]]) == "glogit")
    ##         {
    ##             ## tau <- Mdl.par[[CompCpl]][["tau"]]
    ##             ## tau.a <- 0 ## The lower bound of generalized logit link
    ##             ## tau.b <- 2^(1/2-1/(2*tau)) ## the upper bound

    ##             ## The lower and upper bounds of generalized logit link
    ##             lambdaL <- Mdl.par[[CompCpl]][["lambdaL"]]

    ##             tau.a <- log(2)/(log(2)-log(lambdaL))
    ##             linkCurr$a <- tau.a

    ##             ## tau.b <- linkCurr[["b"]] ## NOTE: Numerical stable. keep it slightly away
    ##             ## from 1.
    ##             ## linkCurr$b <- tau.b

    ##         }
    ##         else
    ##         {
    ##             stop("Such conditional linkage is not implemented!")
    ##         }

    ##         Mdl.par[[CompCpl]][["tau"]] <- parMeanFun(X = XCurr,
    ##                                                 beta = betaCurr,
    ##                                                 linkArgs = linkCurr)
    ##     }
    ## }


    out <- Mdl.par

    return(out)
}
