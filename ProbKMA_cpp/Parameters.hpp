#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <vector>
#include <memory>

enum class stopCriterion{max,mean,quantile};

struct Parameters
{
    Parameters(unsigned int K,const std::vector<unsigned int>& c,
               std::vector<unsigned int>& c_max,unsigned int iter_max,
               stopCriterion stop_criterion,unsigned int quantile,
               double tol,unsigned int iter4elong,double tol4elong,
               double max_elong,unsigned int trials_elong,
               double deltaJK_elong,double max_gap,
               unsigned int iter4clean,double tol4clean,
               double quantile4clean,bool return_options,
               bool parallel);

    unsigned int _K;
    std::vector<unsigned int> _c;
    std::vector<unsigned int> _c_max;
    unsigned int _iter_max;
    stopCriterion _stop_criterion;
    unsigned int _quantile;
    double _tol;
    unsigned int _iter4elong;
    double _tol4elong;
    double _max_elong;
    unsigned int _trials_elong;
    double _deltaJK_elong;
    double _max_gap;
    unsigned int _iter4clean;
    double _tol4clean;
    double _quantile4clean;
    bool _parallel;

}

#endif // PARAMETERS_HPP