logLikelihood <- function(CplNM, Mdl.Y, Mdl.par, Mdl.u, Mdl.d, parUpdate, MCMCUpdateStrategy)
{

###----------------------------------------------------------------------------
### UPDATE MARGINAL PDF AND/OR CDF THE MARGINAL U AND ONLY UPDATED IF THE CORRESPONDING
### PARAMETERS ARE UPDATED.
###----------------------------------------------------------------------------
  CompNM <- names(Mdl.par)
  browser()
  ## Allocate the output structure: Margins + Copula (NA)
  ## Mdl.logLik <- cbind(Mdl.d, NA)
  ## colnames(Mdl.logLik) <- CompNM

  CompUpNM <- unlist(lapply(parUpdate, function(x) any(unlist(x) == TRUE)))
  MargisNM <- CompNM[CompNM  != CplNM]
  MargisUpNM <- CompNM[(CompNM  != CplNM) & CompUpNM]


  ## THE COPULA LOG LIKELIHOOD
  if(tolower(MCMCUpdateStrategy) == "joint")
    {
      evalCpl <- TRUE
      Mdl.PostComp <- lapply(parUpdate, function(x) TRUE)
    }
  else if(tolower(MCMCUpdateStrategy) == "twostage" |
          tolower(MCMCUpdateStrategy) == "margin")
    {
      if(length(MargisUpNM) == 0)
        {
          ## Stage two of the two stage approach
          evalCpl <- TRUE
          Mdl.PostComp <- lapply(parUpdate, function(x) FALSE)
          Mdl.PostComp[[CplNM]] <- TRUE
        }
      else
        {
          evalCpl <- FALSE
          Mdl.PostComp <- lapply(parUpdate, function(x) any(unlist(x) == TRUE))
        }
    }
  else
    {
      stop(paste("MCMC update strategy:", MCMCUpdateStrategy,
                 "not implemented!"))
    }



  ## THE MARGINAL LIKELIHOODS
  densCaller <- list()
  for(CompCaller in MargisNM)
    {
      if(CompCaller %in% MargisUpNM)
        {
          if(tolower(MCMCUpdateStrategy) == "joint")
            {
              densCaller[[CompCaller]] <- c("u", "d")
            }
          else if(tolower(MCMCUpdateStrategy) == "twostage")
            {
              ## Stage two of the two stage approach
              densCaller[[CompCaller]] <- c("u", "d")
            }
          else if(tolower(MCMCUpdateStrategy) == "margin")
            {
              densCaller[[CompCaller]] <- c(NA, "d")
            }
          else
            {
              stop(paste("MCMC update strategy:", MCMCUpdateStrategy,
                         "not implemented!"))
            }
        }
      else
        {
          if(evalCpl == TRUE && any(is.na(Mdl.u)))
            {
              ## if(tolower(MCMCUpdateStrategy) == "joint")
              ##     {
              densCaller[[CompCaller]] <- c("u", "d")
              ##     }
              ## else
              ##     {
              ##         browser()
              ##         stop("Two-stage updating without providing \"u\" is not possible!")
              ##     }
            }
        }
    }

###----------------------------------------------------------------------------
### UPDATING THE LIKELIHOOD
###----------------------------------------------------------------------------

  ## Updating the marginal models TODO: parallel version
  for(CompCaller in names(densCaller))
    {
      Mdl.ud <- MargiModel(
              y = Mdl.Y[[CompCaller]],
              type = MargisTypes[which(MargisNM == CompCaller)],
              par = Mdl.par[[CompCaller]],
              densCaller = densCaller[[CompCaller]])

      if("u" %in% densCaller[[CompCaller]])
        {
          Mdl.u[,  CompCaller] <- Mdl.ud[["u"]]
        }
      if("d" %in% densCaller[[CompCaller]])
        {
          Mdl.d[,  CompCaller] <- Mdl.ud[["d"]]
        }
    }


  ## Updating the copula model if required
  if(evalCpl == TRUE)
    {
      Mdl.logLikCpl <- logCplLik(
              u = Mdl.u,
              CplNM = CplNM,
              parCplRep = Mdl.par[[CplNM]],
              sum = FALSE) # n-by-1

      Mdl.d[, CplNM] <- Mdl.logLikCpl
    }

  out <- list(Mdl.d = Mdl.d, Mdl.u = Mdl.u, Mdl.PostComp = Mdl.PostComp)
  return(out)
}

## This is the optimization version
logLikelihoodOptim <- function(ipar, i, chainCaller, CplNM, Mdl.Y,
                               Mdl.par, Mdl.u, Mdl.d, parUpdate,
                               MCMCUpdateStrategy)
{

  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par[[CompCaller]][[parCaller]][i] <- parVec

  Mdl.ud <- logLikelihood(CplNM = CplNM,
                          Mdl.Y = Mdl.Y,
                          Mdl.par = Mdl.par,
                          Mdl.u = Mdl.u,
                          Mdl.d = Mdl.d,
                          parUpdate = parUpdate,
                          MCMCUpdateStrategy = MCMCUpdateStrategy)
  Mdl.d <- Mdl.ud[["Mdl.d"]]
  Mdl.PostComp <- Mdl.ud[["Mdl.PostComp"]]

  out <- sum(Mdl.d[, unlist(Mdl.PostComp)])
}
