##' Setting up the data object structure.
##'
##' <details>
##' @title <short tile>
##' @param CplNM
##' @param CplParNM
##' @param MargisNM
##' @param MargisParNM
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Tue Jan 03 17:03:24 CET 2012;
##'       Current: Tue Jan 03 17:03:32 CET 2012.
initDataStruc <- function(CplNM, CplParNM, MargisNM, MargisParNM)
  {
    CompNM <- c(MargisNM, CplNM)
    CompParNM <- c(MargisParNM, CplParNM)
    MdlDataStruc <- list()

    for(i in 1:length(CompNM))
        {
            MdlDataStruc[[i]] <- as.list(rep(NA, length(CompParNM[[i]])))
            names(MdlDataStruc[[i]]) <- CompParNM[[i]]
        }
    names(MdlDataStruc) <- CompNM

    out <- MdlDataStruc
    return(out)
  }
