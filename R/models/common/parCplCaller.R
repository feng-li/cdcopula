##' parCplCaller
##'
##' Generate the updating caller.
##' @param parUpdate
##'
##' @param parUpdateOrder
##'
##' @param nIter if nIter is not 1. This indicating the two stage copula estimation is
##' used. Therefore it will generate the complete updating sequence matrix for each
##' iteration.
##'
##' @return
##' @references Li 2014
##' @author Feng Li, Central University of Finance and Economics.
parCplCaller <- function(parUpdate, parUpdateOrder, nIter = 1)
{
    ## check which parameter need to update
    CompNM <- names(parUpdate)
    out.init <- NULL
    for(i in CompNM)
        {
            parNM <- names(parUpdate[[i]])
            parUpdateIdx <- parNM[parUpdate[[i]] == TRUE]
            for(j in parUpdateIdx)
                {

                    out.init <- rbind(out.init , c(i, j))
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
            out <- matrix(out.init, , 2)
        }

    return(out)
}
