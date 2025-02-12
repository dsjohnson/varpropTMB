// namespace varpropTMB {
// 
// template <class Type>
// bool isNA(Type x) {
//   return R_IsNA(asDouble(x));
// }
// 
// template <class Type>
// Type ppois_log(Type x, Type lambda) {
//   return atomic::Rmath::Rf_ppois(asDouble(x), asDouble(lambda), true, true);
// }
// 
// template <class Type>
// Type dstudent(Type x, Type mean, Type sigma, Type df, int give_log = 0) {
//   // from metRology::dt.scaled()
//   // dt((x - mean)/sd, df, ncp = ncp, log = TRUE) - log(sd)
//   Type logres = dt((x - mean) / sigma, df, true) - log(sigma);
//   if (give_log)
//     return logres;
//   else
//     return exp(logres);
// }
// 
// // List of matrices
// template <class Type>
// struct LOM_t : vector<matrix<Type> > {
//   LOM_t(SEXP x) {  // x = list passed from R
//     (*this).resize(LENGTH(x));
//     for (int i = 0; i < LENGTH(x); i++) {
//       SEXP sm = VECTOR_ELT(x, i);
//       (*this)(i) = asMatrix<Type>(sm);
//     }
//   }
// };
// 
// 
// } // namespace end