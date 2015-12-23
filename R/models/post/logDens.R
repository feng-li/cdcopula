logDens <- function(MargisType, Mdl.Y, Mdl.par, Mdl.u, Mdl.d, parUpdate,
                    MCMCUpdateStrategy)
{
###----------------------------------------------------------------------------
### UPDATE MARGINAL PDF AND/OR CDF THE MARGINAL U AND ONLY UPDATED IF THE CORRESPONDING
### PARAMETERS ARE UPDATED.
###----------------------------------------------------------------------------
  if(missing(Mdl.u))
  {
    Mdl.u <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.Y),
                    dimnames = list(NULL, names(Mdl.Y)))
  }

  if(missing(Mdl.d))
  {
    Mdl.d <- matrix(NA, dim(Mdl.Y[[1]])[1], length(Mdl.par),
                    dimnames = list(NULL, names(Mdl.par)))
  }

  CompNM <- names(Mdl.par)
  ## Allocate the output structure: Margins + Copula (NA)
  ## Mdl.logLik <- cbind(Mdl.d, NA)
  ## colnames(Mdl.logLik) <- CompNM

  CplNM <- MargisType[length(MargisType)]
  CompUpNM <- unlist(lapply(parUpdate, function(x) any(unlist(x) == TRUE)))
  MargisNM <- CompNM[CompNM  != CplNM]
  MargisUpNM <- CompNM[(CompNM  != CplNM) & CompUpNM]


  ## THE COPULA LOG LIKELIHOOD
  if(tolower(MCMCUpdateStrategy) == "joint")
  {
    evalCpl <- TRUE
    Mdl.PostComp <- lapply(parUpdate, function(x) TRUE)
  }
  else if(tolower(MCMCUpdateStrategy) == "twostage" ||
          tolower(MCMCUpdateStrategy) == "margin")
  {
    if(length(MargisUpNM) == 0)
    {
      ## Stage two of the two stage approach, only updating copula density
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
        densCaller[[CompCaller]] <- c("u", "d")
      }
    }
  }

###----------------------------------------------------------------------------
### UPDATING THE LIKELIHOOD
###----------------------------------------------------------------------------
  ## Updating the marginal models TODO: parallel version
  for(CompCaller in names(densCaller))
  {
      ## if(CompCaller == "BABA") browser()

    Mdl.ud <- MargiModel(y = Mdl.Y[[CompCaller]],
                         type = MargisType[which(MargisNM == CompCaller)],
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
    Mdl.logLikCpl <- logCplLik(u = Mdl.u,
                               CplNM = CplNM,
                               parCplRep = Mdl.par[[CplNM]],
                               sum = FALSE) # n-by-1
    Mdl.d[, CplNM] <- Mdl.logLikCpl
  }

  out <- list(Mdl.d = Mdl.d, Mdl.u = Mdl.u,
              Mdl.PostComp = Mdl.PostComp)
  return(out)
}
