
#include "ProbKMA.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// definiton of the different distances 
enum class Distance = {L0,L1,L0_L1};

class ProbKMA::_probKMAImp
{
    _probKMAImp(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
                const arma::mat& S0,const Parameters& parameters,
                Distance distance,
                const std::unique_ptr<Dissimilarity>& pdissimilaritY,
                const std::unique_ptr<Motif>& pmotif): _P0(P0), _P1(P1),
                _parameters(parameters)
                {
                    int Y_size = Y.size();
                    int V_size = V.size();
                    if(distance == L0)
                    {   
                        _Y.set_size(Y_size,1);
                        _V.set_size(V_size,1);  
                        for(int i = 0;i<Y_size )
                        {
                            
                        }
                    }
                    if(distance == L2)
                    {

                    }
                    // TODO check lengths are all equal otherwise throw an exception
                    for(int i  = 0; i < enumDistances.size();++i)
                    {
                        _distance = std::move(pdissimilarities[i]);
                        _motif = std::move(pmotifs[i]);
                    
                    }
                }
    
    ~_probKMAImpl() = default;

    // parameters for initialization
    Parameters _parameters;
 
    // different implementations of distance
    std::unique_ptr<Dissimilarity> _distance;
   
    // different implementations of compute_motif
    std::unique_ptr<Motif> _motif;

    Rcpp::List probKMA_run() const
    {
        return Rcpp::List::create();
    }

    void set_parameters(const Parameters& parameters)
    {
        _parameters = paramameters; // check wheter the copy assignment operator is well defined
    }

    void set_distance(const std::unique_ptr<Dissimiilarity>& pnewDistance)
    {
        _distance = std::move(pnewDistance);
    }

    void set_motif(const std::unique_ptr<Motif>& pnewMotif)
    {
        _motif = std::move(pnewMotif);
    }

    // return a R structure
    Rcpp::List toR() const  // da implementare in base ai dati che voglio restituire in R
    {
        return Rcpp::List::create();
    }

    // Functional Data
    arma::field<arma::mat> _Y;
    arma::field<arma::mat> _V;
    arma::mat _P0;
    arma::mat _S0;
    
}


ProbKMA::probKMA(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
                const arma::mat& S0,std::vector<Distance> enumDistances,
                const Parameters& parameters,
                const std::vector<std::unique_ptr<Dissimilarity>>& pdissimilarities,
                const std::vector<std::unique_ptr<Motif>>& pmotif):
                _probKMA(std::make_unique<probKMA_imp>(Y,V,P0,S0,enumDistances,
                                                       parameters,pdissimilarities,pmotif));


Rcpp::List ProbKMA::probKMA_run() const
{
   return _probKMA -> probKMA_run();
}

void ProbKMA::set_parameters(const Parameters& parameters)
{
    _parameters = paramameters;
}

void ProbKMA::set_distance(const std::unique_ptr<Dissimiilarity>& pnewDistance)
{
    _distance = std::move(pnewDistance);
}

void ProbKMA::set_motif(const std::unique_ptr<Motif>& pnewMotif)
{
    _motif = std::move(pnewMotif);
}

