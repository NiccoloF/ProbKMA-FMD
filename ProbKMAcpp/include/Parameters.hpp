#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP
#include <vector>
#include <string>
#include <Rcpp.h>


class Parameters
{
public:
  
    Parameters() = default;
 
    void createParameters(unsigned int K,
                          const std::vector<unsigned int>& c,
                          const std::vector<unsigned int>& c_max,
                          unsigned int iter_max,unsigned int quantile,
                          const std::string& stopCriterion,double tol,
                          unsigned int iter4elong,double tol4elong,
                          double max_elong,unsigned int trials_elong,
                          double deltaJK_elong,double max_gap,
                          unsigned int iter4clean,double tol4clean,
                          double quantile4clean,bool parallel);
    
    void test(SEXP a);

    unsigned int _K;
    std::vector<unsigned int> _c;
    std::vector<unsigned int> _c_max;
    unsigned int _iter_max;
    unsigned int _quantile;
    std::string _stopCriterion;
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
};

#endif // PARAMETERS_HPP