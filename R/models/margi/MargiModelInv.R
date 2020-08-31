#' @export
MargiModelInv <- function(u, par, type)
{
    if(tolower(type) == "gaussian") # continues of course
    {
        mu <- par[["mu"]] # scaler
        phi <- par[["phi"]]   # scaler
        out <- matrix(qnorm(u, mean = mu, sd = phi))
    }
    else if(tolower(type) == "splitt")
    {
        mu <- par[["mu"]] # scaler
        phi <- par[["phi"]]   # scaler
        df <- par[["df"]]   # scaler
        lmd <- par[["lmd"]]   # scaler

        out <- matrix(qsplitt(p = u, mu =  mu, df = df,
                              phi = phi, lmd = lmd))
    }
    else if(tolower(type) %in% c("teigen", "splitt-mixture"))
    {
        ## mu <- par[["mu"]] # scaler
        ## phi <- par[["phi"]]   # scaler
        ## df <- par[["df"]]   # scaler
        ## lmd <- par[["lmd"]]   # scaler
        ## weights = par[["weights"]]

        if(length(par[["lmd"]]) == 0)
        {
            par[["lmd"]] = matrix(1, 1, ncol(par[[1]]))
        }

        out <- matrix(qmixture(p = u, type = "splitt", par.list = par,
                               args = list(nSim = 1000)))

    }
    else if(tolower(type) == "poisson") # the discrete case
    {
        stop("Poisson uses latent variable instead.!")
    }
    else
    {
        stop("The quantile for the marginal density is not implemented yet.")
    }


    ## The output
    return(out)
}
