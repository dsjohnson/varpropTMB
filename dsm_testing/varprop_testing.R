library(Distance)
library(mgcv)
load("~/research/projects/r_packages/varpropTMB/dsm_testing/dsm_test_data.RData")

ddf <- detfc.hr.null
model_data <- dsm.xy.depth$data
formula <- count ~ s(x, y, k=10) + s(depth, k = 20)

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

