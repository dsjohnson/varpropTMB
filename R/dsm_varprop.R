#' @title Fit DSM Model While Accounting for Detection Model Uncertainty
#' @param formula Model formula for the habitat portion of te DSM model.
#' @param seg_data Data frame containing the observations and habitat data for the distance sampling segments
#' @param effort_obj A named list 
#' @param family A GLM or GAM family object
#' 
#' @export

dsm_varprop <- function(formula, model_data, offset_list=NULL, family="poisson", ....){
  # library(mgcv)
  # library(Distance)
  # load("~/research/projects/r_packages/varpropTMB/dsm_testing/dsm_test_data.RData")
  
  ### Design matrices
  sm_list <- parse_smoothers(formula, model_data)
  fe <- all_terms(formula)
  fe <- fe[-get_smooth_terms(fe)]
  if(length(fe)==0) {
    fe_formula <- ~1
  } else{
    fe_formula <- paste0("~ ", paste(fe, collapse=" + "))
  }
  Xf <- model.matrix(fe_formula, model_data)
  y <- model_data[,all.vars(formula)[attr(terms(formula), "response")]]
  
  ### offset
  Oe <- offset_list$offset_func(offset_list$offset_par)
  if(!is.null(offset_list$offset_V)){
    Ve <- offset_list$offset_V
    Ke <- numDeriv::jacobian(offset_list$offset_func, offset_list$offset_par)
    varprop = 1
  } else{
    varprop = 0
    Ve <- 0
    Ke <- 0
  }
  
  ### Lists for dsm model
  dl <- list(
    y=y, Xf=Xf, Xs=do.call(cbind,sm_list$Xs), Zs=sm_list$Zs, 
    Oe=Oe, Ke=Ke, Ve=offset_list$offset_V, 
    family="poisson", n=nrow(model_data), ns=length(sm_list$Zs), 
    varprop = varprop
  )
  
  pl <- list(
    bf = rep(0,ncol(Xf)) , 
    bs = rep(0,ncol(sm_list$Xs)), 
    us = rep(0, sum(sapply(sm_list$Zs, ncol))), 
    ln_sig_us = rep(0,length(sm_list$Zs)), 
    ln_sig_Xs = rep(0,length(sm_list$Xs)), 
    de = ifelse(varprop, rep(0, ncol(Ke)), 0),
    ln_phi = 0, logit_theta = -4.59512
  )
  
  map_list <- lapply(pl, \(x) rep(NA, length(x)))
  # fix family par
  if(family=="poisson"){
    map_list$ln_phi = 0; map_list$logit_theta = -4.59512
  } else if(family=="negbinom"){
    map_list$logit_theta = -4.59512
  }
  # fix smooth
  
  
  
  
} ## end fitting func