library(Distance)
library(mgcv)
library(TMB)
library(sdmTMB)
library(dsm)
library(magrittr)
library(scam)
load("~/research/projects/r_packages/varpropTMB/dsm_testing/dsm_test_data.RData")

ddf <- detfc.hr.null
model_data <- dsm.xy.depth$data
formula <- count ~ t2(x, y, k=10) + s(depth, k = 20, bs="ts")
select <- FALSE
family <- "poisson"

offset_list <- list(
  offset_par = ddf$ddf$par,
  offset_func = function(par){
    ddf$ddf$par <- par
    esw <- predict(ddf, model_data, esw=TRUE, compute = TRUE)$fitted
    off <- log(2 * model_data$Effort * esw)
    return(off)
  },
  offset_V = solve(ddf$ddf$hessian)
)


Ve <- matrix(NA, 3,3)
Ve[2:3,2:3] <- offset_list$offset_V
Ke <- numDeriv::jacobian(offset_list$offset_func, offset_list$offset_par)

model_data$Ke <- Ke


fit <- scam(formula, data=model_data, offset = offset_list$offset_func(offset_list$offset_par),
                 method = "REML", family = quasipoisson(link = "log"))





fit <- sdmTMB(
  count ~ s(x, y, k=10) + s(depth, k = 20),
  data=model_data, family = tweedie(link="log"), spatial = "off",
  offset = log(offset_list$offset_func(offset_list$offset_par)),
  # priors = sdmTMBpriors(
  #   b = mvnormal(location = c(NA,0,0), scale = Ve)
  # )
)

tmb_obj <- fit$tmb_obj
rep <- sdreport(tmb_obj)

preddata$Ke <- matrix(0,nrow(preddata), 2)


preddata$x <- 0
preddata$y <- 0
preddata$depth <- mean(pdata$depth)
ppp <- predict(fit, newdata=preddata, se_fit = TRUE)
ppp$est <- ppp$est-mean(ppp$est)
plot(ppp$depth, ppp$est)

preddata_sf$pred <- ppp$est
mapview(preddata_sf, zcol="pred")

fit_stan <- tmbstan::tmbstan(
  fit$tmb_obj,
  iter = 1500, chains = 1, warmup=500,
  seed = 8675309 # ensures repeatability
)

samps <- sdmTMBextra::extract_mcmc(fit_stan)
pred <- predict(fit, newdata=preddata, mcmc_samples = samps)
preddata_sf$pred <- pred*preddata$area
