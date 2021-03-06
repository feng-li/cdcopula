#' Hamiltonian Monte Carlo with variable selection.
#'
#' Hamiltonian Monte Carlo
#' @param CplNM NA
#' @param Mdl.Y NA
#' @param Mdl.X NA
#' @param Mdl.beta NA
#' @param Mdl.betaIdx NA
#' @param Mdl.parLink NA
#' @param parUpdate NA
#' @param Mdl.priArgs NA
#' @param Mdl.varSelArgs NA
#' @param MCMC.propArgs NA
#' @param Mdl.MargisType NA
#' @param staticCache NA
#' @param Mdl.algorithm NA
#' @param MCMC.UpdateStrategy NA
#' @return "list"
#' @references Li 2012
#' @author Feng Li, Central University of Finance and Economics.
#' @export
HamiltonianMonteCarlo <- function(CplNM, Mdl.Y, Mdl.X, Mdl.beta, Mdl.betaIdx,
                               Mdl.parLink, parUpdate, Mdl.priArgs, Mdl.varSelArgs,
                               MCMC.propArgs, Mdl.MargisType, staticCache, Mdl.algorithm,
                               MCMC.UpdateStrategy)
{
    ## The updating component parameter chain
    chainCaller <- parCplRepCaller(parUpdate)
    CompCaller <- chainCaller[1]
    parCaller <- chainCaller[2]

    ## cat("chainCaller", chainCaller, "\n")

    ## The proposal methods
    algmArgs <- MCMC.propArgs[[CompCaller]][[parCaller]][["algorithm"]]
    beta.MCMC.propArgs <- MCMC.propArgs[[CompCaller]][[parCaller]][["beta"]]
    betaIdx.MCMC.propArgs <- MCMC.propArgs[[CompCaller]][[parCaller]][["indicators"]]

###----------------------------------------------------------------------------
### INITIAL COPY OF CURRENT VALUES
###----------------------------------------------------------------------------

    Mdl.beta.curr <- Mdl.beta
    Mdl.betaIdx.curr <- Mdl.betaIdx
    staticCache.curr <- staticCache

    beta.curr.full <- Mdl.beta[[CompCaller]][[parCaller]] # p-by-lq
    betaIdx.curr <- Mdl.betaIdx[[CompCaller]][[parCaller]] # p-by-lq
    beta.curr <- beta.curr.full[betaIdx.curr] # p*lq

    ## Assume initial no variable selection
    betaIdx.prop <- betaIdx.curr

###----------------------------------------------------------------------------
### VARIABLE SELECTION PROPOSAL
###----------------------------------------------------------------------------

    ## No. of covariates
    nCovs <- dim(betaIdx.curr)[1]
    nPar <- dim(betaIdx.curr)[2]
    ## Randomly propose a subset for covariates to change
    beta01Mat <- matrix(1:(nCovs*nPar), nCovs, nPar)

    varSelCandConfigRow <- Mdl.varSelArgs[[CompCaller]][[parCaller]][["cand"]] # sub.q-by-1

    if(class(varSelCandConfigRow) == "character" &&
       tolower(varSelCandConfigRow) == "2:end")
    {
        varSelCandRow <- 2:nCovs
    }
    else
    {
        varSelCandRow <- varSelCandConfigRow
    }

    if(length(varSelCandRow)>nCovs)
    {
        stop("The number of variable selection candidates should\n",
             "not be greater than the number of covariates for: ",
             paste(chainCaller, collapse = "-"),
             ". Check variable selection settings.")
    }


    varSelCand <- beta01Mat[varSelCandRow, ] # sub.p-by-lq

    ## If variable selection is available, make a change proposal. Otherwise no variable
    ## selection.
    if(length(varSelCand) > 0)
    {
        if(betaIdx.MCMC.propArgs[["type"]] == "binom")
        {
            ## Binomial proposal a small subset
            betaIdx.propCandIdx <- which(rbinom(n = length(varSelCand),
                                                size = 1L,
                                                prob = betaIdx.MCMC.propArgs[["prob"]]) == 1L)

            if(length(betaIdx.propCandIdx) != 0)
            {
                ## Situation when propose a change with prob = betaIdxArgs$prob
                betaIdx.propCand <- varSelCand[betaIdx.propCandIdx]
                betaIdx.prop[betaIdx.propCand] <- !betaIdx.curr[betaIdx.propCand]

                VSProp <- "VSProp"
            }
            else
            {
                ## Situation when propose no changes  =  no variable selection
                VSProp <- NULL
            }
        }
        else
        {
            stop("No such variable selection proposals!")
        }
    }
    else
    {   ## No variable selection
        VSProp <- NULL
    }

    ## The jump density for the variable selection indicators. TODO: Add adaptive scheme

    logJump.Idx.currATprop <- 1
    logJump.Idx.propATcurr <- 1

###----------------------------------------------------------------------------
### INITIALIZING THE MH SCHEMES.
###----------------------------------------------------------------------------

    ## Initially we assume there is only mh step.
    MHUpdate <- "MH"

    ## If variable selection (VS) is proposed, add an extra MH step with VS at the
    ## beginning. otherwise we do the usual update without variable selection,
    ## i.e. betaidx.prop = betaidx.curr.
    MHUpdate <- c(VSProp, MHUpdate)

    ## For each mh step, we obtain an acceptance probability
    nMH <- length(MHUpdate)
    accept.probs <- rep(NA, nMH)

    ## Set a reject flag to handle unexpected situations. If TRUE, the proposal is rejected
    ## anyway regardless of other situations.
    errorFlags <- rep(FALSE, nMH)
###----------------------------------------------------------------------------
### METROPOLIS-HASTINGS ALGORITHM
###----------------------------------------------------------------------------
    for(iMH in 1:nMH)
    {
        if(tolower(algmArgs[["type"]]) == "gnewtonmove")
        { ## Newton method to approach the posterior based on the current draw

            ## cat("beta.curr.full:", beta.curr.full, "\n")

            beta.NTProp <- PropGNewtonMove(Mdl.MargisType = Mdl.MargisType,
                                           MCMC.propArgs = MCMC.propArgs,
                                           Mdl.varSelArgs = Mdl.varSelArgs,
                                           Mdl.priArgs = Mdl.priArgs,
                                           betaIdxProp = betaIdx.prop,
                                           parUpdate = parUpdate,
                                           Mdl.Y = Mdl.Y,
                                           Mdl.X = Mdl.X,
                                           Mdl.parLink = Mdl.parLink,
                                           Mdl.beta = Mdl.beta.curr,
                                           Mdl.betaIdx = Mdl.betaIdx.curr,
                                           staticCache = staticCache.curr,
                                           Mdl.algorithm = Mdl.algorithm,
                                           MCMC.UpdateStrategy = MCMC.UpdateStrategy)

        }
        else if(tolower(algmArgs[["type"]])  == "randomwalk")
        { ## Random walk metropolis (with/without variable selection)
            stop("Not implement yet!")

            ## beta.NTProp <- list(errorFlag = FALSE,
            ##                     param = ??,
            ##                     HessObsInv = ??)
        }
        else
        {
            stop("Not implement yet!")
        }

        ## Check if it is a good proposal
        if(beta.NTProp$errorFlag)
        { ## Something is wrong. Skip this MH update(if is VS, jump to pure update without
            ## variable selection) or terminate the algorithm.

            errorFlags[iMH] <- TRUE
            betaIdx.prop <- betaIdx.curr
            next
        }

        ## else Continues with the Metropolis-Hastings Propose a draw from multivariate
        ## t-distribution based on the proposed An idea (out of loud) : Generate a matrix of
        ## param. Select the one that give max acceptance probability (but need correct the
        ## acceptance probability.)

        ## The information for proposed density via K-step Newton's method
        beta.prop.mean <- matrix(beta.NTProp[["param"]], 1) # 1-by-p
        beta.prop.sigma <- -beta.NTProp[["HessObsInv"]]


        staticCache.prop <- beta.NTProp[["staticCache"]]
        if(tolower(beta.MCMC.propArgs[["type"]]) == "mvt")
        {
            ## require("mvtnorm")
            ## The proposal parameters block
            beta.prop <- (beta.prop.mean +
                          rmvt(sigma = (beta.prop.sigma+t(beta.prop.sigma))/ 2,
                               n = 1, df = beta.MCMC.propArgs[["df"]],
                               method = "chol", checkSymmetry = FALSE))
        }

###----------------------------------------------------------------------------
### REVERSE PROBABILITY OF THE PROPOSALS
###----------------------------------------------------------------------------

        ## Newton method to approach the posterior for the proposed draw
        beta.prop.full <- matrix(0, nCovs, nPar)
        beta.prop.full[betaIdx.prop] <- beta.prop
        Mdl.beta.prop <- Mdl.beta.curr
        Mdl.beta.prop[[CompCaller]][[parCaller]] <- beta.prop.full

        Mdl.betaIdx.prop <- Mdl.betaIdx.curr
        Mdl.betaIdx.prop[[CompCaller]][[parCaller]] <- betaIdx.prop

        if(tolower(algmArgs[["type"]]) == "gnewtonmove")
        {
            ## cat("beta.prop.full:", beta.prop.full, "\n")

            beta.NTPropRev <- PropGNewtonMove(Mdl.MargisType = Mdl.MargisType,
                                              MCMC.propArgs = MCMC.propArgs,
                                              Mdl.varSelArgs = Mdl.varSelArgs,
                                              Mdl.priArgs = Mdl.priArgs,
                                              betaIdxProp = betaIdx.curr,
                                              parUpdate = parUpdate,
                                              Mdl.Y = Mdl.Y,
                                              Mdl.X = Mdl.X,
                                              Mdl.parLink = Mdl.parLink,
                                              Mdl.beta = Mdl.beta.prop,
                                              Mdl.betaIdx = Mdl.betaIdx.prop,
                                              Mdl.algorithm = Mdl.algorithm,
                                              staticCache = staticCache,
                                              MCMC.UpdateStrategy = MCMC.UpdateStrategy)

        }
        else if(tolower(algmArgs[["type"]])  == "randomwalk")
        {
            stop("Not implement yet!")

            ## beta.NTPropRev <- list(errorFlag = FALSE,
            ##                        param = ??,
            ##                        HessObsInv = ??)
        }
        else
        {
            stop("Not implement yet!")
        }

        if(beta.NTPropRev$errorFlag)
        { ## Something is wrong.
            errorFlags[iMH] <- TRUE
            betaIdx.prop <- betaIdx.curr
            next
        }

        ## The information for proposed density via K-step Newton's method
        beta.propRev.mean <- matrix(beta.NTPropRev[["param"]], 1) # 1-by-p
        beta.propRev.sigma <- -beta.NTPropRev[["HessObsInv"]] # p-by-p

        ## browser()
###----------------------------------------------------------------------------
### COMPUTING THE METROPOLIS-HASTINGS RATIO
###----------------------------------------------------------------------------

        ## The jump density for proposed point at proposed mode and the jump density for
        ## current draw at reverse proposed mode.
        if(tolower(beta.MCMC.propArgs[["type"]]) == "mvt")
        {

            logJump.propATprop <- dmvt(x = beta.prop - beta.prop.mean,
                                       sigma = (beta.prop.sigma+t(beta.prop.sigma))/2,
                                       df = beta.MCMC.propArgs[["df"]], log = TRUE,
                                       checkSymmetry = FALSE)

            logJump.currATpropRev<- dmvt(x = beta.curr - beta.propRev.mean,
                                         sigma = (beta.propRev.sigma+t(beta.propRev.sigma))/2,
                                         df = beta.MCMC.propArgs[["df"]], log = TRUE,
                                         checkSymmetry = FALSE)
        }

        ## The log posterior for the proposed draw
        logPost.propOut <- logPost(Mdl.MargisType = Mdl.MargisType,
                                   Mdl.Y = Mdl.Y,
                                   Mdl.X = Mdl.X,
                                   Mdl.beta = Mdl.beta.prop,
                                   Mdl.betaIdx = Mdl.betaIdx.prop,
                                   Mdl.parLink = Mdl.parLink,
                                   Mdl.varSelArgs = Mdl.varSelArgs,
                                   Mdl.priArgs = Mdl.priArgs,
                                   parUpdate = parUpdate,
                                   Mdl.algorithm = Mdl.algorithm,
                                   staticCache = staticCache,
                                   MCMC.UpdateStrategy = MCMC.UpdateStrategy)

        logPost.prop <- logPost.propOut[["Mdl.logPost"]]
        staticCache.prop <- logPost.propOut[["staticCache"]]

        ## The log posterior for the current draw
        logPost.curr <- logPost(Mdl.MargisType = Mdl.MargisType,
                                Mdl.Y = Mdl.Y,
                                Mdl.X = Mdl.X,
                                Mdl.beta = Mdl.beta.curr,
                                Mdl.betaIdx = Mdl.betaIdx.curr,
                                Mdl.parLink = Mdl.parLink,
                                Mdl.varSelArgs = Mdl.varSelArgs,
                                Mdl.priArgs = Mdl.priArgs,
                                Mdl.algorithm = Mdl.algorithm,
                                parUpdate = parUpdate,
                                staticCache = staticCache,
                                MCMC.UpdateStrategy = MCMC.UpdateStrategy)[["Mdl.logPost"]]

        ## Compute the (log) MH ratio.
        logMHRatio <- (logPost.prop - logPost.curr +
                       logJump.currATpropRev - logJump.propATprop +
                       logJump.Idx.currATprop - logJump.Idx.propATcurr)

        ## The acceptance probability

        if(is.na(logMHRatio))
        {
            ## bad proposal, i.e logJump.currATpropRev = -Inf, or logJump.propATprop =
            ## -Inf
            accept.prob.curr <- 0
        }
        else
        {
            accept.prob.curr <- exp(min(0, logMHRatio))
        }

        ## cat(iMH, logPost.prop, logPost.curr, logJump.propATprop, logJump.currATpropRev, accept.prob.curr, "\n")
###----------------------------------------------------------------------------
### THE MH ACCEPTANCE PROBABILITY AND ACCEPT/REJECT THE PROPOSED DRAW.
###----------------------------------------------------------------------------

        if(runif(1) < accept.prob.curr) #!is.na(accept.prob.curr)
        {
            ## keep the proposal
            betaIdx.curr <- betaIdx.prop
            beta.curr <- beta.prop

            Mdl.beta.curr <- Mdl.beta.prop
            Mdl.betaIdx.curr <- Mdl.betaIdx.prop

            staticCache.curr <- staticCache.prop
            ## browser()
        }
        else
        {## keep the current
            betaIdx.prop <- betaIdx.curr
        }

        accept.probs[iMH] <- accept.prob.curr
    }

###----------------------------------------------------------------------------
### THE FINAL OUTPUT
###----------------------------------------------------------------------------
    ## The acceptance prob are from the last MH update or keep current draw.
    if(errorFlags[nMH] == TRUE)
    { ## epic failure
        out <- list(errorFlag = TRUE)
    }
    else
    {
        out <- list(beta =  Mdl.beta.curr[[CompCaller]][[parCaller]],
                    betaIdx = Mdl.betaIdx.curr[[CompCaller]][[parCaller]],
                    accept.prob = accept.probs[nMH],
                    staticCache = staticCache.curr,
                    errorFlag = FALSE)
    }
    return(out)
}


#' Leapfrog function in Hamiltonian Monte Carlo
#'
#'
#' @param par The parameter vector
#' @param phi The momentum vector which has the same length of the parameter vector.
#' @param Mass_inv The inverse of mass matrix, the covariance of the momentum distribution
#'     which is usually set as a diagonal.
#' @param epsilon Scaler, Scale factor.
#' @param gradfun The gradient function for the log posterior with respect to its
#'     parameter vector.
#' @param ... Other parameters passed to `gradfun`.
#' @return list
#' @references Bayesian Data Analysis , 3rd,  p 301. Hoffman & Gelman (2014)
#' @author Feng Li, Central University of Finance and Economics.
#' @export
LeapFrog <- function(par, momentum, Mass_inv, epsilon, gradfun, ...)
{
    momentum_new <- momentum + 0.5 * epsilon * gradfun(par, ...)
    par_new <- par + epsilon * Mass_inv %*% momentum_new
    momentum_new2 <- momentum_new + 0.5 * epsilon * gradfun(par_new, ...)
    out <- list(par = par_new, momentum = momentum_new2)
    return(out)
}


#' Heuristic for choosing an initial value of epsilon in HMC
#'
#' @param par NA
#' @param f NA
#' @param gradfun NA
#' @param Mass NA
#' @param epsilon NA
#' @param maxIter NA
#' @return NA
#' @references NA
#' @author Feng Li, Central University of Finance and Economics.
#' @export
FindEpsilon = function(par, f, gradfun, Mass, epsilon = 1, maxIter = 100)
{
    momentum <- rmvnorm(length(par), 0, Mass)

    Mass_inv <- ginv(Mass)
    proposed <- LeapFrog(par = par,
                         momentum = momentum,
                         Mass_inv = Mass_inv,
                         epsilon = epsilon,
                         gradfun = gradfun, ...)
    log_ratio <- (joint_log_density(proposed$par, proposed$r, f, Mass_inv) -
                  joint_log_density(par, r, f, Mass_inv))

    alpha <- ifelse(exp(log_ratio) > 0.5, 1, -1)
    if(is.nan(alpha)) alpha <- -1

    iIter <- 1
    while(is.nan(log_ratio) || iIter <= maxIter ||
          (alpha * log_ratio > (-alpha)*log(2)))
    {
        epsilon <- 2^alpha * epsilon
        proposed <- LeapFrog(par = par,
                             momentum = momentum,
                             Mass_inv = Mass_inv,
                             epsilon = epsilon,
                             gradfun = gradfun, ...)
        log_ratio <- (joint_log_density(proposed$par, proposed$r, f, Mass_inv) -
                      joint_log_density(par, r, f, Mass_inv))

        iIter <- iIter + 1
        if(iIter  ==  maxIter)
        {
            message("Could not find reasonable epsilon in ", maxIter, " iterations!")
        }
    }
    message("Epsilon = ", epsilon, " found after ", iIter, " steps.")
    return(epsilon)
}
