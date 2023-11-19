#include "Parameters.hpp"

Parameters::Parameters(unsigned int K,const std::vector<unsigned int>& c,
                      std::vector<unsigned int>& c_max,unsigned int iter_max,
                      stopCriterion stop_criterion,unsigned int quantile,
                      double tol,unsigned int iter4elong,double tol4elong,
                      double max_elong,unsigned int trials_elong,
                      double deltaJK_elong,double max_gap,
                      unsigned int iter4clean,double tol4clean,
                      double quantile4clean,bool return_options,
                      bool parallel):_K(K),_c(c),_c_max(std::move(c_max)),
                      _iter_max(iter_max),_stop_criterion(stop_criterion),
                      _quantile(quantile),_tol(tol),_iter4elong(iter4elong),
                      _tol4elong(tol4elong),_max_elong(max_elong),
                      _trials_elong(trials_elong),_deltaJK_elong(deltaJK_elong),
                      _max_gap(max_gap),_iter4clean(iter4clean),
                      _tol4clean(tol4clean),_quantile4clean(quantile4clean),
                      _parallel(parallel) {};
