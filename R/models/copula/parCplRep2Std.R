##' <description>
##'
##' <details>
##' @title <short tile>
##' @param copula
##' @param parCplRep
##' @return
##' @references
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: ; Current: .
parCplRep2Std <- function(CplNM, parCplRep)
  {
      ## Initialize output Storage and name it
      out <- list()

      if(tolower(CplNM) == "bb7")
          {
              ## The reparameterized parameters input
              lambdaL <- parCplRep[["lambdaL"]]
              tau <- parCplRep[["tau"]]

              lambdaU <- as.vector(kendalltauInv(
                  CplNM = CplNM, parCplRep = parCplRep))

              ## The standard parametrization
              delta <- -log(2)/log(lambdaL)
              theta <- log(2)/log(2-lambdaU)

              ## The first parameter
              out[["delta"]] <- delta

              ## The second parameter
              out[["theta"]] <- theta
          }
      else if(tolower(CplNM) == "mvt")
          {
              lambdaL <- parCplRep[["lambdaL"]]
              tau <- parCplRep[["tau"]]

              rho <- sin(tau*pi/2) # n-by-lq

              FUN <- function(x1, x2)
                  {
                      parCpl <- list(rho = x1, df = x2)
                      lambda(CplNM = "mvt", parCpl = parCpl)[["lambdaL"]]
                  }

              df <- funinv2d(FUN = FUN,
                             method = "tabular",
                             x1 = rho,
                             y = lambdaL,
                             x1lim = c(0, 1),
                             x2lim = c(2, 30),
                             tol = 1e-3) # n-by-1

              out[["rho"]] <- rho
              out[["df"]] <- df
          }
    else
      {
        stop("This copula is not implemented!")
      }
    return(out)
  }
