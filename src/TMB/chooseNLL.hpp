/// @file chooseNLL.hpp

#ifndef chooseNLL_hpp
#define chooseNLL_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

enum valid_fam {
  gaus_fam = 1,
  gamm_fam = 2
};

template<class Type>
Type nll(Type dat, Type a, Type b, int fam)
{
  Type out;
  switch(fam) {
  case gaus_fam:
    out = -dnorm(dat,a,b,true);
    break;
  case gamm_fam:
    out =  -dgamma(dat,a,b,true);
    break;
  default:
    error("Nope!");
  }
  return out;
}


/// Negative log-likelihood of the chosen distribution.
template<class Type>
Type chooseNLL(objective_function<Type>* obj) {
  DATA_INTEGER(n);
  DATA_VECTOR(x); // data vector
  DATA_INTEGER(fam); //family
  PARAMETER_VECTOR(par); // mean parameter
  Type out = 0;
  Type a = par(0);
  Type b = par(1);
  for(int i=0; i<n; i++) {
    out += nll(x(i), a, b, fam);
  }
  return out; // negative log likelihood
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
