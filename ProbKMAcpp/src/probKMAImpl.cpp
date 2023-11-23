#include "ProbKMA.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

class ProbKMA::_probKMAImp
{
public:
    _probKMAImp(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
                const arma::mat& S0,SEXP parameters,
                SEXP dissimilarity,SEXP motif): _P0(P0), _S0(S0)
                {
                    set_parameters(parameters);
                    set_distance(dissimilarity);
                    set_motif(motif);
                    Initialize(Y,V);
                }
    
    ~_probKMAImp() = default;

    // parameters for initialization
    Parameters* _parameters;
 
    // different implementations of distance
    Dissimilarity* _distance;
   
    // different implementations of compute_motif
    Motif* _motif;

    void Initialize(const Rcpp::List& Y,const Rcpp::List& V)
    {
        int Y_size = Y.size();
        int V_size = V.size();
        if(!Rf_isNull(Y[0]) and !Rf_isNull(Y[1]))
        {
            const Rcpp::List& Y0_temp = Y[0];
            const Rcpp::List& Y1_temp = Y[1];
            const Rcpp::List& V0_temp = V[0];
            const Rcpp::List& V1_temp = V[1];
            _Y.set_size(Y0_temp.size(),2);
            _V.set_size(V0_temp.size(),2);
            for(int i = 0; i < Y0_temp.size(); i++)
            {
                _Y(i,0) = Rcpp::as<arma::mat>(Y0_temp[i]);
                _Y(i,1) = Rcpp::as<arma::mat>(Y1_temp[i]);
            }
            for(int i = 0; i < V0_temp.size(); i++)
            {
                _V(i,0) = Rcpp::as<arma::mat>(V0_temp[i]);
                _V(i,1) = Rcpp::as<arma::mat>(V1_temp[i]); 
            }
        }
        else if(!Rf_isNull(Y[0]))
        {
            const Rcpp::List& Y0_temp = Y[0];
            const Rcpp::List& V0_temp = V[0];
            _Y.set_size(Y0_temp.size(),1);
            _V.set_size(V0_temp.size(),1);
            for(int i = 0; i < Y0_temp.size(); i++)
            {
                _Y(i,0) = Rcpp::as<arma::mat>(Y0_temp[i]);
            }
            for(int i = 0; i < V0_temp.size(); i++)
            {
                _V(i,0) = Rcpp::as<arma::mat>(V0_temp[i]);
            }
        }
        else
        {
            const Rcpp::List& Y1_temp = Y[1];
            const Rcpp::List& V1_temp = V[1];
            _Y.set_size(Y1_temp.size(),1);
            _V.set_size(V1_temp.size(),1);
            for(int i = 0; i < Y1_temp.size(); i++)
            {
                _Y(i,0) = Rcpp::as<arma::mat>(Y1_temp[i]);
            }
            for(int i = 0; i < V1_temp.size(); i++)
            {
                _V(i,0) = Rcpp::as<arma::mat>(V1_temp[i]);
            }
        }
    } 

    Rcpp::List probKMA_run() const
    {
      
        return Rcpp::List::create();
    }

    void set_parameters(SEXP newParameters)
    {
      Rcpp::S4 parametersObj(newParameters);
      Rcpp::Environment env(parametersObj);
      Rcpp::XPtr<Parameters> xptr( env.get(".pointer") );
      _parameters = static_cast<Parameters*> (R_ExternalPtrAddr(xptr));
    }
    
    void set_distance(SEXP newDistance)
    {
      Rcpp::S4 distanceObj(newDistance);
      Rcpp::Environment env(distanceObj);
      Rcpp::XPtr<Dissimilarity> xptr( env.get(".pointer") );
      _distance = static_cast<Dissimilarity*> (R_ExternalPtrAddr(xptr));
    }

    void set_motif(SEXP newMotif)
    {
      Rcpp::S4 motifObj(newMotif);
      Rcpp::Environment env(motifObj);
      Rcpp::XPtr<Motif*> xptr( env.get(".pointer") );
      _motif = static_cast<Motif*> (R_ExternalPtrAddr(xptr));
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
  
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,const arma::mat& P0,
                 const arma::mat& S0,SEXP parameters,
                 SEXP dissimilarity,SEXP motif):
                 _probKMA(std::make_unique<_probKMAImp>(Y,V,P0,S0,parameters,
                                                       dissimilarity,motif)) {}; 


Rcpp::List ProbKMA::probKMA_run() const
{
   return _probKMA -> probKMA_run();
}

void ProbKMA::set_parameters(SEXP newParameters)
{
    _probKMA -> set_parameters(newParameters);
}

void ProbKMA::set_distance(SEXP newDistance)
{
  _probKMA -> set_distance(newDistance);
}

void ProbKMA::set_motif(SEXP newMotif)
{
    _probKMA -> set_motif(newMotif);
}

