#' parCplRepCaller
#'
#' Generate the updating caller.
#' @param parUpdate NA
#' @param parUpdateOrder NA
#' @param MCMC.nIter if MCMC.nIter is not 1. This indicating the two stage copula estimation is
#' used. Therefore it will generate the complete updating sequence matrix for each
#' iteration.
#'
#' @return NA
#' @references Li 2014
#' @author Feng Li, Central University of Finance and Economics.
#' @export
parCplRepCaller <- function(parUpdate, parUpdateOrder)
{
    CompNM <- names(parUpdate)

    ## check which parameter need to update
    out.init <- NULL
    for(CompCaller in CompNM)
    {
        parNM <- names(parUpdate[[CompCaller]])
        parUpdateIdx <- parNM[parUpdate[[CompCaller]] == TRUE]
        for(parCaller in parUpdateIdx)
        {
            out.init <- rbind(out.init , c(CompCaller, parCaller))
        }
    }

    ## Order the parameters if asked
    if(!missing(parUpdateOrder))
    {
        parOrder <- unlist(parUpdateOrder)
        parUpdate <- unlist(parUpdate)
        parUpdateOrderVec <- order(parOrder[parUpdate])
        out <- out.init[parUpdateOrderVec, , drop = FALSE]
    }
    else
    {
        out <- matrix(out.init, length(out.init)/2, 2)
    }

    return(out)
}


#' @export
parCplRepCallerVine <- function(parUpdate)
{
    ## check which parameter need to update
    out.init <- NULL
    for(CompCaller in 1:length(parUpdate))
    {
        ## parNM <- names(parUpdate[[CompCaller]])
        parLst <- 1:length(parUpdate[[CompCaller]])

        parUpdateIdx <- parLst[parUpdate[[CompCaller]] == TRUE]
        for(parCaller in parUpdateIdx)
        {
            out.init <- rbind(out.init , c(CompCaller, parCaller))
        }
    }

    out <- matrix(out.init, length(out.init)/2, 2)

    return(out)
}
