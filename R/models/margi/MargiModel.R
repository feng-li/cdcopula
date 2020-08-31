#' CDF and PDF of the marginal distribution.
#'
#' The detailed description can be found in the main setting file for each
#' input variable.
#' @param y "vector".
#'        The response variables for the marginal model.
#' @param type "vector" with character input.
#'        The type of the marginal model.
#' @param par "list".
#'        The parameters input for the marginal model.
#' @return "list".
#'
#'        Return the percentile and log density of the marginal models.
#'
#'
#' @references Li 2012
#' @author Feng Li, Department of Statistics, Stockholm University, Sweden.
#' @note Created: Tue Jan 17 19:27:25 CET 2012; Current: Mon Jan 05 16:24:51 CST 2015.
#' @export
MargiModel <- function(y, type, par, densCaller = c("u", "d"))
{
###----------------------------------------------------------------------------
### CONVERT FOREIGN MODELS INTO STANDARD SPECIFICATIONS
###----------------------------------------------------------------------------
    if(all((names(par) %in% c("mu", "phi"))) && length(names(par)) == 2)
    {
        typeStd <- "gaussian"
    }
    else if(all((names(par) %in% c("mu", "phi", "df", "lmd"))) &&
            length(names(par)) == 4)
    {
        typeStd <- "splitt"
    }
    else
    {
        typeStd <- type
    }


###----------------------------------------------------------------------------
### THE STANDARD MARGINAL MODELS
###----------------------------------------------------------------------------
    if(tolower(typeStd) %in% c("gaussian"))
    {
        ## The mean and standard deviation for Gaussian density
        mu <- par[["mu"]] # mean parameter
        phi <- par[["phi"]]   # standard deviation

        ## The percentile representation
        if("u" %in% tolower(densCaller))
        {
            u <- pnorm(y, mean = mu, sd = phi, log = FALSE)
        }
        ## The quantile representation
        if("d" %in% tolower(densCaller))
        {
            d <- dnorm(y, mean = mu, sd = phi, log = TRUE)
        }
    }
    else if (tolower(typeStd)  == "splitt")
    {
        ## The marginal likelihood
        ## Literal translation from GSMMatlab code AsymStudT

        mu <- par[["mu"]]  # location parameter
        df <- par[["df"]] # Degrees of freedom
        phi <- par[["phi"]] # Scaling Parameter
        lmd <- par[["lmd"]] # Skewness Parameter

        ## CDF
        if("u" %in% tolower(densCaller))
        {
            u <- psplitt(q = y, mu = mu, df = df, phi = phi, lmd = lmd)
        }
        ## PDF
        if("d" %in% tolower(densCaller))
        {
            d <- dsplitt(x = y, mu = mu, df = df, phi = phi, lmd = lmd,
                         log = TRUE)
        }
    }
    else if (tolower(typeStd)  == "teigen")
    {
        ## The marginal likelihood
        ## Literal translation from GSMMatlab code AsymStudT
        par[["lmd"]] <- array(1, dim(par[["weights"]])) # not skewed t

        ## CDF
        if("u" %in% tolower(densCaller))
        {
            u <- pmixture(q = y, type = "splitt", par.list = par,  log = FALSE)
        }
        ## PDF
        if("d" %in% tolower(densCaller))
        {
            d <- dmixture(x = y, type = "splitt", par.list = par,  log = TRUE)
        }

    }
    else if (tolower(typeStd)  == "poisson")
    {
        mu <- par[["mu"]]  # location parameter

        if("u" %in% tolower(densCaller))
        {
            ## Reserve space
            u <- matrix(0, length(y), 2, dimnames = list(NULL, c("u", "um")))

            ## Stober, et al CSDA 2015. No need `du` or `dm`
            up <- ppois(q = y, lambda = mu, log.p = FALSE)

            Idx <- (y>=1) # valid idx to avoid lower boundary for calculating `um`

            ym <- y[Idx]-1

            um <- ppois(q = ym, lambda = mu[Idx], log.p = FALSE)

            ## Assign results
            u[, 1] <- up
            u[Idx, 2] <- um
        }

        if("d" %in% tolower(densCaller))
        {
            d <- dpois(x = y, lambda = mu, log = TRUE)
        }

    }
    else
    {
        stop("This type of margin is not implemented.")
    }


    ## The out storage
    out <- list()

    ## The output
    if("u" %in% tolower(densCaller))
    {
        ## Numeric stability check for u

        tol <- 1e-6
        u0.idx <- (u<tol)
        if(any(u0.idx))
        {
            u[u0.idx] <- 0+tol
            warning(paste("Marginal CDF too close to 0, replaced with", tol))
        }

        u1.idx <- ((1-u)<tol)
        if(any(u1.idx))
        {
            u[u1.idx] <- 1-tol
            warning(paste("Marginal CDF too close to 1, replaced with", 1-tol))
        }

        ## assign column names
        if(ncol(u) ==  1)
        {
            colNM <- "u"
        }
        else
        {
            colNM <- c("u", "um")
        }
        colnames(u) <- colNM


        out[["u"]] <- u

    }
    if("d" %in% tolower(densCaller))
    {
        out[["d"]] <- d

    }
    return(out)
}
