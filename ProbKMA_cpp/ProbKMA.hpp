#ifndef PROBKMA_HPP
#define PROBKMA_HPP

#include "RcppArmadillo.h"
#include "Parameters.hpp"
#include "Dissimilarity.hpp"
#include "Motif.hpp"
#include <vector>
#include <ranges>
#include <algorithm>
#include <memory>
#include <unordered_map>

// Create a hash function for the enum
struct EnumClassHash {
    template <typename T>
    std::size_t operator()(T t) const {
        return static_cast<std::size_t>(t);
    }
};

// Forward declaration
enum class Distance: unsigned int;
class _probKMAImp;


class ProbKMA
{
 public:
    probKMA(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
            const arma::mat& S0,const Parameters& parameters,
            Distance distance,
            const std::unique_ptr<Dissimilarity>& pdissimilaritY,
            const std::unique_ptr<Motif>& pmotif);
    
    ~probKMA() = default;
    
    // run probKMA's algorithm
    Rcpp::List probKMA_run() const;

    void set_distance(const std::unique_ptr<Dissimiilarity>& pnewDistance);
    void set_motif(const std::unique_ptr<Motif>& pnewMotif);

 private:

    // Pimpl design
    class _probKMAImp;
    std::unique_ptr<probKMA_imp> _probKMA;
}

#endif // PROBKMA_HPP




