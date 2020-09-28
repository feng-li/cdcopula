#' @export
logLikSubSampling <- function(dens,Mdl.Y,Mdl.beta, nObsRaw)
{
    browser()
    nObsUsed <-length(Mdl.Y[[1]])

    ##calculate the mean and variance of sample##
    mu <- lapply(Mdl.Y, mean)
    sigma <-lapply(Mdl.Y, var)

###control variate qi####
#####Assuming a norm function,try other forms#####
    mean_SML<-mean(Mdl.beta[[1]][1][[1]])
    mean_OEX<-mean(Mdl.beta[[2]][1][[1]])
    q1<-dnorm(mu[1],mean_SML,sd=sigma,log=TRUE)
    l1<-q1
############part1   l_m##########
#####depend on qi and after-transfered dens_i###

    l_m<-(dens[1]-q1)/Mdl.dataUsedIdx/nObsRaw+q1

    for (i in 2:50)
    {
        q<-dnorm(mu[i],mean_SML,sd=sigma,log=TRUE)
        l<-(dens[i]-q)/n/nObsRaw
        l_m<-l_m+q+l
    }

########part2   sigma_m###############
#######depend on qi ,li and l_m#######
    sigma_m<-((l1-q1)/Mdl.dataUsedIdx-l_m)^2

    for (i in 2:Mdl.dataUsedIdx)
    {
        q <-dnorm(mu[i],mean_SML,sigma,log=TRUE)
        l<-(dens[i]-q)/n
        sigma_m<-sigma_m+(l-l_m)^2
    }
    sigma_m<-sigma_m/nObsRaw/(nObsRaw-1)

######combine#######
    logLikSumSampling<-exp(l_m-sigma_m^2/2)

    ## the final output###
    out<-logLikSumSampling
}
