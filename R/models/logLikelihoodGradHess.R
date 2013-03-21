##' The gradient for the log posterior in the copula model
##'
##' This function usually does not need to modify. Note that this function does
##' not give the complete posteriors but its components,  i.e. gradient for the
##' likelihood function, gradient for the prior,  gradient for the linkage.
##' @param CplNM "character".
##' @param MargisTypes "list".
##' @param Mdl.Y "list".
##' @param Mdl.X "list".
##' @param Mdl.parLink "list".
##' @param Mdl.beta "list".
##' @param Mdl.betaIdx "list".
##' @param parUpdate "list".
##' @param varSelArgs "list".
##' @param staticCache "list".
##' @param gradMethods "character"
##'
##'        If is "numeric", the numeric gradient is returned; else return the
##' analytical gradient.
##'
##' @return "list".
##'
##'        errorFlag: If anything went wrong, quit with TRUE.
##'
##'        logLilGradObs: Gradient for the log likelihood
##'
##' @references Li 2012
##' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
##' @note Created: Thu Feb 02 22:45:42 CET 2012;
##'       Current: Thu Feb 02 22:45:48 CET 2012.
logLikelihoodGradHess <- function(
    CplNM,
    MargisTypes,
    Mdl.Y,
    Mdl.X,
    Mdl.parLink,
    Mdl.beta,
    Mdl.betaIdx,
    parUpdate,
    varSelArgs,
    staticCache,
    gradMethods = c("analytic", "numeric")[1])
{


  ## The updating chain
  chainCaller <- parCplCaller(parUpdate)
  CompCaller <- chainCaller[1]
  parCaller <- chainCaller[2]

  Mdl.par <- staticCache[["Mdl.par"]]
  ## if(is(Mdl.par, "try-error")) browser()

  ## if(parCaller == "tau") browser()

###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE MARGINAL LIKELIHOOD
###----------------------------------------------------------------------------

  if(tolower(CompCaller) != tolower(CplNM))
    {
      yCurr <- Mdl.Y[[CompCaller]]
      parCurr <- Mdl.par[[CompCaller]]
      typeCurr <- MargisTypes[[CompCaller]]

      ## Gradient Fraction in the marginal component. n-by-1
      if("analytic" %in% tolower(gradMethods))
        {
          MargiGradObs.ana <- MargiModelGrad(
              par = parCurr,
              y = yCurr,
              type = typeCurr,
              parCaller = parCaller)
          MargiGradObs <- MargiGradObs.ana
        }
      if("numeric" %in% tolower(gradMethods))
        {
          ## A simple wrapper that calculates the numerical gradient for given
          ## parameters

          ## NOTE: that this part depends on Chris Sims' numgrad() function. It
          ## is under the flutils module.

          MargiModelGradNumFun <- function(
              x, parCaller, parCurr, yCurr, typeCurr)
            {
              parCurr[[parCaller]] <- x
              MargiLogLikObs <- MargiModel(
                  par = parCurr,
                  y = yCurr,
                  type = typeCurr)$u
              out <- MargiLogLikObs
              return(out)
            }

          nObs <- length(Mdl.Y[[1]])
          MargiGradObs.num <- matrix(NA, nObs, 1)
          for(i in 1:nObs)
            {
              gradTry <- try(numgrad(
                  fcn = MargiModelGradNumFun,
                  x = parCurr[[parCaller]][i],
                  parCaller = parCaller,
                  parCurr = lapply(parCurr, function(x, i)x[i], i = i),
                  yCurr = yCurr[i],
                  typeCurr = typeCurr), silent = TRUE)

              if(is(gradTry, "try-error"))
                {
                  MargiGradObs.num[i] <- NA
                }
              else
                {
                  MargiGradObs.num[i] <- gradTry$g
                }
            }
          MargiGradObs <- MargiGradObs.num
        }

      ## plot(MargiGradObs.ana, MargiGradObs.num)

      ## The Copula parameter caller is the marginal CDF, i.e. u1,  u2, ...
      cplCaller <- paste("u", which(names(MargisTypes)%in%CompCaller), sep = "")

      staticCache[["Mdl.u"]][, CompCaller] <- MargiModel(
          y = yCurr,
          type = typeCurr,
          par = parCurr)[["u"]]
    }
  else
    {
      ## Only update the gradient for copula parameters
      ## Gradient Fraction in the copula component.
      MargiGradObs <- 1
      cplCaller <- parCaller
    }

  ## Error checking
  if(any(is.na(MargiGradObs)) || any(is.infinite(MargiGradObs)))
    {
      return(list(errorFlag = TRUE))
    }


###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE COPULA LIKELIHOOD
###----------------------------------------------------------------------------


  if("analytic" %in% tolower(gradMethods))
    {

      ## The gradient for the copula function. n-by-1
      logCplGradObs.ana <- logCplGrad(
          CplNM = CplNM,
          u = staticCache$Mdl.u,
          parCpl = Mdl.par[[CplNM]],
          cplCaller = cplCaller,
          Mdl.X = Mdl.X,
          Mdl.beta = Mdl.beta)
      logCplGradObs <- logCplGradObs.ana
    }
  if("numeric" %in% tolower(gradMethods))
    {
      ## The gradient for the copula function. scaler input and output
      ## NOTE: The numerical gradient may not work well if the tabular version
      ## of Kendall's tau is used due to the precision
      logCplGradNumFun <- function(x, u,  CompCaller, parCaller, cplCaller,
                                   CplNM, parCpl, staticCache)
        {
          if(tolower(cplCaller) %in% c("u1", "u2"))
            {
              ## Calling the marginal CDF u1, u2
              ## u <- staticCache$Mdl.u
              u[, CompCaller] <- x
            }
          else
            {
              ## Calling copula parameters
              parCpl[[parCaller]] <- x
            }
          out <- logCplLik(u = u, CplNM = CplNM,
                           parCpl = parCpl,
                           logLik = FALSE)
          return(out)
        }

      nObs <- length(Mdl.Y[[1]])
      logCplGradObs.num <- matrix(NA, nObs, 1)
      for(i in 1:nObs)
        {
          if(tolower(cplCaller) %in% c("u1", "u2"))
            {
              ## Calling the marginal CDF u1, u2
              xCurr <- staticCache$Mdl.u[i, CompCaller]
            }
          else
            {
              ## Calling copula parameters
              xCurr <- Mdl.par[[CompCaller]][[parCaller]][i]
            }

          ## if(i == 144 && parCaller == "lambdaL") browser()

          gradTry <- try(numgrad(
              fcn = logCplGradNumFun,
              x = xCurr,
              u = staticCache$Mdl.u[i, , drop = FALSE],
              CompCaller = CompCaller,
              parCaller = parCaller,
              cplCaller = cplCaller,
              CplNM =  CplNM,
              parCpl = lapply(Mdl.par[[CplNM]], function(x, i)x[i], i = i),
              staticCache = staticCache), silent = TRUE)
          if(is(gradTry, "try-error"))
            {
              logCplGradObs.num[i] <- NA
            }
          else
            {
              logCplGradObs.num[i] <- gradTry$g
            }

        }
      ## The Gradient For The Link Function n-by-1
      logCplGradObs <- logCplGradObs.num
    }

  if(parCaller %in% c("tau", "lambdaL"))
    {
      ## browser()
      ## cat(parCaller, Mdl.par[[CompCaller]][[parCaller]][1], "\n")
      ## plot(logCplGradObs.ana, logCplGradObs.num)
    }
  ## Error checking
  if(any(is.na(logCplGradObs)) || any(is.infinite(logCplGradObs)))
    {
      return(list(errorFlag = TRUE))
    }



###----------------------------------------------------------------------------
### GRADIENT FRACTION IN THE LINK FUNCTION
###----------------------------------------------------------------------------
  ## The gradient for the link function n-by-1
  LinkGradObs <- parCplMeanFunGrad(
      CplNM = CplNM,
      Mdl.par = Mdl.par,
      Mdl.parLink = Mdl.parLink,
      chainCaller = chainCaller)

  ## Error checking
  if(any(is.na(LinkGradObs)) || any(is.infinite(LinkGradObs)))
    {
      return(list(errorFlag = TRUE))
    }

###----------------------------------------------------------------------------
### THE OUTPUT
###----------------------------------------------------------------------------

  ## The gradient for the likelihood,  n-by-1
  logLikGradObs <- (logCplGradObs*MargiGradObs)*LinkGradObs

  ## par(mfrow = c(1, 2))

  ## The output
  out <- list(logLikGradObs = logLikGradObs,
              logLikHessObs = NA,
              errorFlag = FALSE)


  plot(logLikGradObs)
  print(chainCaller)
  if(chainCaller[2] == "df") browser()

  return(out)
}
