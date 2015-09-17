


require("rmgarch")
spec <- gogarchspec(mean.model  =  list(model  =  "constant",  robust  = FALSE,
                                        lag  =  1,  lag.max  =  NULL,
                                        lag.criterion  = "AIC",
                                        external.regressors  =  NULL,
                                        robust.control  =  list("gamma"  =  0.25,
                                                                "delta"  =  0.01,
                                                                "nc" =  10,
                                                                "ns"  =  500)),

                    variance.model  =  list(model  =  "sGARCH",
                                            garchOrder  =  c(1, 1),
                                            submodel  =  NULL,
                                            variance.targeting  =  FALSE),
                    distribution.model  =  "mvnorm")

data <- cbind(Y[[1]], Y[[2]])

fit <- gogarchfit(spec = spec,  data = data, out.sample = 0)
resid <- as.ts(fit@mfit$residuals)


pred <- gogarchforecast(fit = fit, n.ahead = 10)


plot(resid)
