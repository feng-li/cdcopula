lambdaL.mean <- read.table("~/Desktop/lambdaLmean.txt", header = TRUE)
library(tseries)# package for unitroots test
colnames(lambdaL.mean)[-1] = rep(c('None', 'Macro', 'Specific', 'Macro+Specific'), 3)
SZKFT.CGWC = ts(lambdaL.mean[, 2:5],start = c(2005, 1), frequency = 4)

SZKFT.CGWC.df=fortify(SZKFT.CGWC)
colnames(SZKFT.CGWC.df)[1] = 'Time'
mdf <- melt(SZKFT.CGWC.df, id.vars="Time", value.name="value", variable.name="(SZKFT, CGWC)")
p1 <- ggplot(data=mdf, aes(x=Time, y=value, group = `(SZKFT, CGWC)`, colour = `(SZKFT, CGWC)`, linetype=`(SZKFT, CGWC)`)) +
  geom_line(size = 0.8) + scale_linetype_manual(values=c("solid", "dotted", "dashed","dotdash"))+
  ylim(0, 1) + labs(y = 'Tail dependence', x = '') + 
  theme(text = element_text(size=15)) +
  scale_x_date(breaks = date_breaks("1 year"), date_labels = '%Y') +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))

SZKFT.XITG = ts(lambdaL.mean[, 6:9],start = c(2005, 1), frequency = 4)
SZKFT.XITG.df=fortify(SZKFT.XITG)
colnames(SZKFT.XITG.df)[1] = 'Time'
mdf <- melt(SZKFT.XITG.df, id.vars="Time", value.name="value", variable.name="(SZKFT, XITG)")
p2 <- ggplot(data=mdf, aes(x=Time, y=value, group = `(SZKFT, XITG)`, colour = `(SZKFT, XITG)`, linetype=`(SZKFT, XITG)`)) +
  geom_line(size = 0.8) + scale_linetype_manual(values=c("solid", "dotted", "dashed","dotdash"))+
  ylim(0, 1) + labs(y = 'Tail dependence', x = '') + 
  theme(text = element_text(size=15))+
  scale_x_date(breaks = date_breaks("1 year"), date_labels = '%Y')+
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))

CGWC.XITG = ts(lambdaL.mean[, 10:13],start = c(2005, 1), frequency = 4)
CGWC.XITG.df=fortify(CGWC.XITG)
colnames(CGWC.XITG.df)[1] = 'Time'
mdf <- melt(CGWC.XITG.df, id.vars="Time", value.name="value", variable.name="(CGWC, XITG)")
p3 <- ggplot(data=mdf, aes(x=Time, y=value, group = `(CGWC, XITG)`, colour = `(CGWC, XITG)`, linetype=`(CGWC, XITG)`)) +
  geom_line(size = 0.8) + scale_linetype_manual(values=c("solid", "dotted", "dashed","dotdash"))+
  labs(y = 'Tail dependence') + 
  theme(text = element_text(size=15))+
  scale_x_date(breaks = date_breaks("1 year"), date_labels = '%Y')+
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1))

library(gridExtra)
grid.arrange(p1, p2, p3, ncol = 1)
