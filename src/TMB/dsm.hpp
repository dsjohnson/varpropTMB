/// @file dsm.hpp

#ifndef dsm_hpp
#define dsm_hpp

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR obj

// negative log-likelihood of the gamma distribution
template<class Type>
Type dsm(objective_function<Type>* obj) {
  DATA_VECTOR(x); // data vector
  PARAMETER(alpha); // shape parameter
  PARAMETER(beta); // scale parameter
  return -sum(dgamma(x, alpha, beta, true));
}

#undef TMB_OBJECTIVE_PTR
#define TMB_OBJECTIVE_PTR this

#endif
