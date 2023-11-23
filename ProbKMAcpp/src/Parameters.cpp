#include "Parameters.hpp"

void Parameters::createParameters(unsigned int K,
                                  const std::vector<unsigned int>& c,
                                  const std::vector<unsigned int>& c_max,
                                  unsigned int iter_max,unsigned int quantile,
                                  const std::string& stopCriterion,double tol,
                                  unsigned int iter4elong,double tol4elong,
                                  double max_elong,unsigned int trials_elong,
                                  double deltaJK_elong,double max_gap,
                                  unsigned int iter4clean,double tol4clean,
                                  double quantile4clean,bool parallel)
  {
    _K = K;
    _c = c;
    _c_max = c_max;
    _iter_max = iter_max;
    _quantile = quantile;
    _stopCriterion = stopCriterion;
    _tol = tol;
    _iter4elong = iter4elong;
    _tol4elong = tol4elong;
    _max_elong = max_elong;
    _trials_elong = trials_elong;
    _deltaJK_elong = deltaJK_elong;
    _max_gap = max_gap;
    _iter4clean = iter4clean;
    _tol4clean = tol4clean;
    _quantile4clean = quantile4clean;
    _parallel = parallel;
  }

  // Expose the custom constructor to R
  RCPP_MODULE(ParametersModule) {
    Rcpp::class_<Parameters>("Parameters")
    .constructor()
    .method("build", &Parameters::createParameters)
    .method("test", &Parameters::test)
    .field("K",&Parameters::_K)
    .field("c",&Parameters::_c)
    .field("c_max",&Parameters::_c_max)
    .field("iter_max",&Parameters::_iter_max)
    .field("quantile",&Parameters::_quantile)
    .field("stopCriterion",&Parameters::_stopCriterion)
    .field("tol",&Parameters::_tol)
    .field("iter4elong",&Parameters::_iter4elong)
    .field("tol4elong",&Parameters::_tol4elong)
    .field("max_elong",&Parameters::_max_elong)
    .field("trials_elong",&Parameters::_trials_elong)
    .field("deltaJK_elong",&Parameters::_deltaJK_elong)
    .field("max_gap",&Parameters::_max_gap)
    .field("iter4clean",&Parameters::_iter4clean)
    .field("tol4clean",&Parameters::_tol4clean)
    .field("quantile4clean",&Parameters::_quantile4clean)
    .field("parallel",&Parameters::_parallel);
  }

