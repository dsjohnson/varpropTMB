/// @file dsm.hpp

using namespace density;

#ifndef dsm_hpp
#define dsm_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// Family spec
// enum valid_family {
//   pois_family = 1,
//   nb_family = 2,
//   tw_family = 3
// };



// negative log-likelihood of the gamma distribution
template<class Type>
Type dsm(objective_function<Type>* obj) {
  
  // Data objects
  
  // DATA_INTEGER(family);
  DATA_INTEGER(n);
  DATA_VECTOR(y);
  DATA_MATRIX(Xf);
  // DATA_INTEGER(null_sel);
  DATA_VECTOR(Oe);
  
  // DATA_INTEGER(has_smooth);
  // DATA_INTEGER(n_bs);
  DATA_INTEGER(n_us);
  DATA_MATRIX(Xs);
  // DATA_IVECTOR(Xs_idx);
  DATA_MATRIX(Zs);
  DATA_IVECTOR(Zs_idx);
  
  // DATA_INTEGER(varprop);
  DATA_MATRIX(Ke);
  DATA_MATRIX(Ve);
  
  // Parameters
  PARAMETER_VECTOR(bf);
  PARAMETER_VECTOR(bs);
  PARAMETER_VECTOR(us);
  PARAMETER_VECTOR(ln_sig_us);
  // PARAMETER_VECTOR(ln_sig_bs);
  PARAMETER_VECTOR(de);
  // PARAMETER(ln_phi);
  // PARAMETER(logit_p);
  

  Type jnll = 0.0;
  
  // Transform parameters
  vector<Type> sig_us = exp(ln_sig_us);
  // vector<Type> sig_bs = exp(ln_sig_bs);
  // Type p = invlogit(logit_p) + 1;
  // Type phi = exp(ln_phi);
  
  // Add smooths
  //if(has_smooth) 
  // ln_mu += Xs*bs + Zs*us;
  
  // Varprop
  //if(varprop) 
  //ln_mu += Ke*de;
    
  vector<Type> ln_mu = Oe + Xf*bf +  Xs*bs + Zs*us + Ke*de;
  
  
  // Inverse link (has to be log) of linear predictor
  vector<Type> mu = exp(ln_mu);
  
  
  // obs_likelihood
  for(int i=0; i<n; i++){
  //   switch(family) {
  //   case pois_family:
      jnll -= dpois(y(i), mu(i), true);
  //     break;
  //   case nb_family: {
  //     Type s1 = ln_mu(i); // log(mu_i)
  //     Type s2 = 2. * s1 - ln_phi; // log(var - mu)
  //     jnll -= dnbinom_robust(y(i), s1, s2, true);
  //   }
  //     break;
  //   case tw_family:
  //     jnll -= dtweedie(y(i), mu(i), phi, p, true);
  //     break;
  //   default:
  //     error("Unrecognized family!");
  //   }
  }
  
  
  // Add smooth penalty
  // if(has_smooth){
    for(int j=0; j<n_us; j++) jnll -= dnorm(us(j), Type(0), sig_us(Zs_idx(j)), true);
  // }
  // if(has_smooth * null_sel){
    //for(int j=0; j<n_bs; j++) jnll -= dnorm(bs(j), Type(0), sig_bs(Xs_idx(j)), true);
  // }
  
  // Add varprop penalty
  // if(varprop) 
    jnll += MVNORM(Ve)(de);
  
  
  return jnll;
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
