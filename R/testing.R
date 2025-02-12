
#' @title choose nll
#' @importFrom TMB MakeADFun
#' @export
choose_nll <- function(x, a, b, fam=1){
  data_list <- list(model="chooseNLL", n=length(x), x=x, fam=fam)
  param_list <- list(par=c(a,b))
  ll <- TMB::MakeADFun(data = data_list,
                       parameters = param_list,
                       DLL = "varpropTMB_TMBExports")
  
  return(list(
    ll$fn(param_list$par), # call its negative loglikelihood
    ll$gr(param_list$par)
  )
  )
}