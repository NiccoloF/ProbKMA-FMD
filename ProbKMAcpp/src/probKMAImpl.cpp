#include "ProbKMA.hpp"

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

class ProbKMA::_probKMAImp
{
public:
  
    _probKMAImp(const Rcpp::List& Y,const Rcpp::List& V,const Rcpp::List& params,
                const arma::mat& P0,const arma::mat& S0)
                :_parameters(params),_P0(P0),_S0(S0)
                {
                    Initialize(Y,V);
                }

   
    ~_probKMAImp() = default;

    // parameters for initialization
    
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

    Rcpp::List probKMA_run(const SEXP& dissimilarity,const SEXP& motif) const
    {
      Rcpp::S4 dissimilairtyObj(dissimilarity);
      Rcpp::S4 motifObj(motif);
      Rcpp::Environment env_diss(dissimilairtyObj);
      Rcpp::Environment env_motif(motifObj);
      Rcpp::XPtr<Dissimilarity> xptr_diss( env_diss.get(".pointer") );
      Rcpp::XPtr<MotifBase> xptr_motif( env_motif.get(".pointer") );
      Dissimilarity* diss_ptr = static_cast<Dissimilarity*> (R_ExternalPtrAddr(xptr_diss));
      MotifBase* motif_ptr = static_cast<MotifBase*> (R_ExternalPtrAddr(xptr_motif));
      
      
        return Rcpp::List::create();
    }

    void set_parameters(const Rcpp::List& newParameters)
    {
      _parameters = Parameters(newParameters);
    }
    
    
    // return a R structure
    Rcpp::List toR() const  // da implementare in base ai dati che voglio restituire in R
    {
        return Rcpp::List::create();
    }
  
    // Functional Data
    Parameters _parameters;
    arma::field<arma::mat> _Y;
    arma::field<arma::mat> _V;
    arma::mat _P0;
    arma::mat _S0;
  
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,
                 const Rcpp::List& parameters,
                 const arma::mat& P0,const arma::mat& S0):
                 _probKMA(std::make_unique<_probKMAImp>(Y,V,parameters,P0,S0)) {}; 


Rcpp::List ProbKMA::probKMA_run(const SEXP& dissimilarity,const SEXP& motif) const
{
   return _probKMA -> probKMA_run(dissimilarity,motif);
}

void ProbKMA::set_parameters(const Rcpp::List& newParameters)
{
    _probKMA -> set_parameters(newParameters);
}


RCPP_EXPOSED_CLASS(ProbKMA);

RCPP_MODULE(ProbKMAModule) {
  Rcpp::class_<ProbKMA>("ProbKMA")
  .constructor<Rcpp::List,Rcpp::List,
               Rcpp::List,arma::mat,
              arma::mat>()
  .method("run", &ProbKMA::probKMA_run)
  .method("set_parameters", &ProbKMA::set_parameters);
}
