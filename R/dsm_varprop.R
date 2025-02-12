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
  source("~/research/projects/r_packages/varpropTMB/R/smooth_processing.R", echo=FALSE)
  
  ### Design matrices
  # Response
  y <- model_data[,all.vars(formula)[attr(terms(formula), "response")]]
  
  # Fixed effects
  fe <- all_terms(formula)
  fe <- fe[-get_smooth_terms(fe)]
  if(length(fe)==0) {
    fe_formula <- ~1
  } else{
    fe_formula <- paste0("~ ", paste(fe, collapse=" + "))
  }
  Xf <- model.matrix(fe_formula, model_data)
  
  ### Smooths
  source("~/research/projects/r_packages/varpropTMB/R/smooth_processing.R", echo=FALSE)
  sm_list <- parse_smoothers(formula, model_data)
  if(sm_list$has_smooths){
    Xs_idx <- NULL
    for(i in 1:length(sm_list$Xs)){
      Xs_idx <- c(Xs_idx, rep(i-1, ncol(sm_list$Xs[[i]])))
    }
    Xs <- do.call(cbind,sm_list$Xs)
    Zs_idx <- NULL
    for(i in 1:length(sm_list$Zs)){
      Zs_idx <- c(Zs_idx, rep(i-1, ncol(sm_list$Zs[[i]])))
    }
    Zs <- do.call(cbind,sm_list$Zs)
  } else {
    Xs_idx <- Xs <- Zs_idx <- Zs <- 0
  }
 
  ### Effort offset
  Oe <- offset_list$offset_func(offset_list$offset_par)
  if(!is.null(offset_list$offset_V)){
    Ve <- offset_list$offset_V
    Ke <- numDeriv::jacobian(offset_list$offset_func, offset_list$offset_par)
    varprop <- 1
  } else{
    varprop <- 0
    Ve <- 0
    Ke <- 0
  }
  
  family_int <- which(family==c("poisson","nbinom","tweedie"))
  
  ### Lists for dsm model
  
  dl <- list(
    model="dsm", family = family_int, y=y, n=length(y), Xf=Xf, null_sel=as.integer(select),
    Oe=Oe, has_smooth=as.integer(sm_list$has_smooths), n_bs=length(Xs_idx),
    n_us=length(Zs_idx), Xs=Xs, Xs_idx=Xs_idx, Zs=Zs, Zs_idx=Zs_idx,
    varprop=varprop, Ke=Ke, Ve=offset_list$offset_V, 
    family="poisson"
  )
  
  pl <- list(
    bf = rep(0,ncol(dl$Xf)) , 
    bs = rep(0,ncol(dl$Xs)), 
    us = rep(0, length(dl$Zs_idx)), 
    ln_sig_us = rep(0,max(dl$Zs_idx)+1), 
    ln_sig_bs = rep(0,max(dl$Xs_idx)+1), 
    de = ifelse(varprop, rep(0, ncol(dl$Ke)), 0),
    ln_phi = 0, logit_p = -4.59512
  )
  
  ml <- lapply(pl, \(x) seq_along(x))
  # fix family par
  if(family=="poisson"){
    ml$ln_phi <- NA; ml$logit_p<- NA
  } else if(family=="nbinom"){
    ml$logit_p <- NA
  }
  # fix smooth
  if(!dl$has_smooth){
    ml$bs <- NA; ml$us <- NA; ml$ln_sig_us <- NA; ml$ln_sig_bs <- NA
  }
  # fix select
  if(!select & sm_list$has_smooths){
    ml$ln_sig_bs <- rep(NA,max(Xs_idx)+1)
  }
  
  ml <- lapply(ml, factor)
  
  # dyn.load(dynlib("/Users/devin.johnson/research/projects/r_packages/varpropTMB/src/TMB/varpropTMB_TMBExports"))
  obj <- TMB::MakeADFun(data = dl, parameters = pl, map = ml, random = c("de","us"),
                        DLL = "varpropTMB_TMBExports")
  
  opt <- nlminb(obj$par, obj$fn, gradient = obj$gr, control=list(iter.max=10000, eval.max=10000))
  
  obj <- TMB::MakeADFun(data = dl, parameters = pl, map = ml, 
                        DLL = "varpropTMB_TMBExports")
  
  
  
} ## end fitting func