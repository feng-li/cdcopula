setwd("~/code/cdcopula/inst/extdata/CreditRisk/CEI/Y")

DTD <- read.table(file = "DD_data.txt", header = TRUE)
colnames(DTD) <- c("date", "quarter", "SZKFT", "SZSED", "CGWC", "HDEIT", "GWII", "SHBL", "CNSS", "XITG", "CECC", "NJPE")

library("GGally")

order <- c(3, 5, 10, 4, 6, 7, 8, 9, 11, 12)

DTDOrdered = ts(DTD[, order], start = c(2005, 1), frequency = 4)
autoplot(DTDOrdered, ylab = 'Distance to Default') +
  annotate("text", x = 2012, y = -2.2, label = "CGWC ", size = 3, col = 'black')+
  annotate("text", x = 2013.1, y = -0.5, label =  " XITG", size = 3, col = 'black')+
  annotate("text", x = 2011.9, y = 0, label = "SZKFT ", size = 3, col = 'black') +
  theme(legend.title = element_blank()) +
  scale_y_continuous(breaks=seq(-2.5, 0.8, 0.5)) +
  scale_x_continuous(breaks=seq(2005, 2016, 1))

ggpairs(DTD[, order],
        upper  =  list(continuous  =  wrap("density", alpha  =  0.5)),
        lower  =  list(continuous  =  wrap("points",  size = 0.8))) +
        theme_grey(base_size  =  8) +
        theme(axis.text.x = element_text(angle = 90))



rm(list = ls())
load("~/running/CreditRisk/R/000021_000066_macroeconomic_covariate_BB7_copula_model.Rdata")
sourceDir("~/code/cdcopula/R", recursive = TRUE)

OUT.TS <- CplMCMC.summary(OUT.FITTED[[1]])

lambdaL.mean <- OUT.TS[["par.summary"]][["ts.mean"]][["BB7"]][["lambdaL"]]
lambdaL.sd <- OUT.TS[["par.summary"]][["ts.sd"]][["BB7"]][["lambdaL"]]


rm(list = ls())
setwd("~/code/cdcopula/inst/figs/credit-risk")
lambdaL.mean <- read.table("lambdaLmean.txt", header = TRUE)
library(tseries)# package for unitroots test
SZKFT.CGWC.lambdaLmean.None <- ts(lambdaL.mean[, 2], start = c(2005, 1), frequency = 4)
SZKFT.CGWC.lambdaLmean.Macro <- ts(lambdaL.mean[, 3], start = c(2005, 1), frequency = 4)
SZKFT.CGWC.lambdaLmean.Spe <- ts(lambdaL.mean[, 4], start = c(2005, 1), frequency = 4)
SZKFT.CGWC.lambdaLmean.MacroSpe <- ts(lambdaL.mean[, 5], start = c(2005, 1), frequency = 4)

SZKFT.XITG.lambdaLmean.None <- ts(lambdaL.mean[, 6], start = c(2005, 1), frequency = 4)
SZKFT.XITG.lambdaLmean.Macro <- ts(lambdaL.mean[, 7], start = c(2005, 1), frequency = 4)
SZKFT.XITG.lambdaLmean.Spe <- ts(lambdaL.mean[, 8], start = c(2005, 1), frequency = 4)
SZKFT.XITG.lambdaLmean.MacroSpe <- ts(lambdaL.mean[, 9], start = c(2005, 1), frequency = 4)

CGWC.XITG.lambdaLmean.None <- ts(lambdaL.mean[, 10], start = c(2005, 1), frequency = 4)
CGWC.XITG.lambdaLmean.Macro <- ts(lambdaL.mean[, 11], start = c(2005, 1), frequency = 4)
CGWC.XITG.lambdaLmean.Spe <- ts(lambdaL.mean[, 12], start = c(2005, 1), frequency = 4)
CGWC.XITG.lambdaLmean.MacroSpe <- ts(lambdaL.mean[, 13], start = c(2005, 1), frequency = 4)

#pdf("tail_dependence.pdf", family = "GB1")
par(mfrow = c(3, 1))
plot(SZKFT.CGWC.lambdaLmean.None,ylim=c(0,1),col="black", lty = "dashed")
lines(SZKFT.CGWC.lambdaLmean.Macro, col = "red", lty = "dotted")
lines(SZKFT.CGWC.lambdaLmean.Spe, col = "purple", lty = "dash-dotted")
lines(SZKFT.CGWC.lambdaLmean.MacroSpe, col = "green")

plot(SZKFT.XITG.lambdaLmean.None, ylim = c(0,1), col = "black")
lines(SZKFT.XITG.lambdaLmean.Macro, col = "red")
lines(SZKFT.XITG.lambdaLmean.Spe, col = "purple")
lines(SZKFT.XITG.lambdaLmean.MacroSpe, col = "green")

plot(CGWC.XITG.lambdaLmean.None, ylim = c(0,1), col = "black", ylab = "lambdaL", xlab = "time")
lines(CGWC.XITG.lambdaLmean.Macro, col = "red", ylab = "lambdaL", xlab = "time")
lines(CGWC.XITG.lambdaLmean.Spe, col = "purple", ylab = "lambdaL", xlab = "time")
lines(CGWC.XITG.lambdaLmean.MacroSpe, col = "green", ylab = "lambdaL", xlab = "time")
