#include "ProbKMA.hpp"
#include "Factory.hpp"
#include <forward_list>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]


class ProbKMA::_probKMAImp
{
public:
  
    _probKMAImp(const Rcpp::List& Y,const Rcpp::List& V,
                const Rcpp::List& parameters,
                const matrix& P0,const imatrix& S0,
                const std::string& diss):
                _parameters(parameters),_P0(P0),_S0(S0)
               {  
                    // Initialize c++ Data structure
                    Initialize(Y,V);
                    
                    // Create Dissimilarity factory
                    util::SharedFactory<Dissimilarity> dissfac;
                    dissfac.FactoryRegister<L2>("L2",_parameters._w);
                    dissfac.FactoryRegister<H1>("H1",_parameters._w,_parameters._alpha);
                    
                    //Create Motif factory
                    util::SharedFactory<MotifPure> motfac;
                    motfac.FactoryRegister<Motif_L2>("L2");
                    motfac.FactoryRegister<Motif_H1>("H1");
                    
                    //check whether diss is valid
                    _motfac = motfac.instantiate(diss); // copy elision
                    _dissfac = dissfac.instantiate(diss);
                    if(not(_motfac and _dissfac))
                      Rcpp::Rcerr<<"Invalid dissimilarity: Choose between L2,H1"<<std::endl;
                }

    ~_probKMAImp() = default;

    void Initialize(const Rcpp::List& Y,const Rcpp::List& V)
    {
        // Convert Rcpp Data Structure(List) into Armadillo data Structure(field)
        if(!Rf_isNull(Y[0]) and !Rf_isNull(Y[1]))
        {
            const Rcpp::List& Y0 = Y[0];
            const Rcpp::List& Y1 = Y[1];
            const Rcpp::List& V0 = V[0];
            const Rcpp::List& V1 = V[1];
            std::size_t Y_size = Y0.size();
            std::size_t V_size = V0.size();
            _Y.set_size(Y_size,2);
            _V.set_size(V_size,2);
            for(int i = 0; i < Y_size; i++)
            {
                _Y(i,0) = Rcpp::as<matrix>(Y0[i]);
                _Y(i,1) = Rcpp::as<matrix>(Y1[i]);
            }
            for(int i = 0; i < V_size; i++)
            {
                _V(i,0) = Rcpp::as<matrix>(V0[i]);
                _V(i,1) = Rcpp::as<matrix>(V1[i]); 
            }
        }
        else if(!Rf_isNull(Y[0]))
        {
            const Rcpp::List& Y0 = Y[0];
            const Rcpp::List& V0 = V[0];
            std::size_t Y_size = Y0.size();
            std::size_t V_size = V0.size();
            _Y.set_size(Y_size,1);
            _V.set_size(V_size,1);
            for(int i = 0; i < Y_size; i++)
            {
                _Y(i,0) = Rcpp::as<matrix>(Y0[i]);
            }
            for(int i = 0; i < V_size; i++)
            {
                _V(i,0) = Rcpp::as<matrix>(V0[i]);
            }
        }
        else
        {
            const Rcpp::List& Y1 = Y[1];
            const Rcpp::List& V1 = V[1];
            std::size_t Y_size = Y1.size();
            std::size_t V_size = V1.size();
            _Y.set_size(Y_size,1);
            _V.set_size(V_size,1);
            for(int i = 0; i < Y_size; i++)
            {
                _Y(i,0) = Rcpp::as<matrix>(Y1[i]);
            }
            for(int i = 0; i < V_size; i++)
            {
                _V(i,0) = Rcpp::as<matrix>(V1[i]);
            }
        }
    } 

    Rcpp::List probKMA_run() 
    {
      /// Iterate ////////////////////////////////////
      std::size_t iter = 0;
      std::forward_list<double> J_iter;
      std::forward_list<double> BC_dist_iter;
      auto BC_dist = std::numeric_limits<double>::infinity();
      const unsigned int& iter_max = _parameters._iter_max;
      const vector quantile4clean(_parameters._quantile4clean); // Necessary to compute quantile
      matrix D;
      umatrix keep;
      Rcpp::Rcout<<"Hello0"<<std::endl;
      std::size_t _n_rows_V = _V.n_rows;
      std::vector<uvector> V_dom(_n_rows_V);
      arma::field<arma::mat> V_new(_n_rows_V,2);
      Rcpp::Rcout<<"Hello1"<<std::endl;
      while(iter < iter_max and BC_dist > _parameters._tol)
      {
        iter++;
        ///// clean motifs //////////////////////////
        const auto P_old = _P0;
        const unsigned int& iter4clean = _parameters._iter4clean;
        const double& tol4clean = _parameters._tol4clean;
        if((iter>1)&&(!(iter%iter4clean))&&(BC_dist<tol4clean))
        {
          keep = D < arma::quantile(D,quantile4clean);
          const uvector empty_k = arma::find(arma::sum(keep,0));
          if(!empty_k.empty())
          {
            arma::uvec index_min = arma::index_min(D,0);
            keep(empty_k(index_min),index_min).fill(1);
          }
          _P0.zeros();
          _P0.elem(arma::find(keep)).fill(1); // set one where values of keep are true
        }
        Rcpp::Rcout<<"Hello2"<<std::endl;
        for(int i = 0;i < _n_rows_V;++i)
        {
          Rcpp::Rcout<<"Hello10"<<std::endl;
          Rcpp::Rcout<<"V_size = "<<_V(i,0).size()<<std::endl;
          const uvector& V_dom_temp = util::findDomain<matrix>(_V(i,0));
          // capire se ha senso dichiararlo fuori
          Rcpp::Rcout<<"Hello11"<<std::endl;
          Rcpp::Rcout<<"V_dom_temp_size"<<V_dom_temp.size()<<std::endl;
          const auto V_new_variant = _motfac->compute_motif(V_dom_temp,_S0.col(i),
                                                            _P0.col(i),_Y,
                                                            _parameters._m);
          
          Rcpp::Rcout<<"Hello12"<<std::endl;
          if(auto ptr_1 = std::get_if<MotifPure::indexField>(&V_new_variant))
          {
            const arma::sword& index = ptr_1->second;
            _S0.col(i) += index;
            V_new.row(i) = ptr_1->first; 
          }
          else
          {
            Rcpp::Rcout<<"Hello4"<<std::endl;
            V_new.row(i) = *(std::get_if<arma::field<arma::mat>>(&V_new_variant));
          }
          V_dom[i] = util::findDomain<matrix>(V_new(i,0));
        }
        
      }
      Rcpp::Rcout<<"Hello5"<<std::endl;
      return Rcpp::List::create(Rcpp::Named("P0")=_P0,Rcpp::Named("S0")=_S0,
                                Rcpp::Named("V_dom")=V_dom);
    }
  
    void set_parameters(const Rcpp::List& newParameters)
    {
      _parameters = newParameters;
    }
    
    
    // return a R structure
    Rcpp::List toR() const  // da implementare in base ai dati che voglio restituire in R
    {
        return Rcpp::List::create();
    }
  
    // Functional Data
    arma::field<matrix> _Y;
    arma::field<matrix> _V;
    
    //Motif and dissimilarity
    std::unique_ptr<MotifPure> _motfac;
    std::unique_ptr<Dissimilarity> _dissfac;
    
    //Parameters
    Parameters _parameters;
    
    // Membership and shifting matrix
    matrix _P0;
    imatrix _S0;
    
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,
                 const Rcpp::List& parameters,
                 const matrix& P0,const imatrix& S0,
                 const std::string& diss):
                 _probKMA(std::make_unique<_probKMAImp>(Y,V,parameters,P0,S0,diss)) {}; 


Rcpp::List ProbKMA::probKMA_run() const
{
   return _probKMA -> probKMA_run();
}

void ProbKMA::set_parameters(const Rcpp::List& newParameters)
{
    _probKMA -> set_parameters(newParameters);
}

// [[Rcpp::export(initialChecks)]]
Rcpp::List initialChecks(const Rcpp::List& Y0,const Rcpp::List& Y1,
                         const Rcpp::NumericMatrix& P0,
                         const Rcpp::NumericMatrix& S0,
                         const Rcpp::List& params,
                         const Rcpp::String& diss,
                         const double alpha,
                         const Rcpp::NumericVector& w)
{
  try 
  {
    Rcpp::Environment base("package:ProbKMAcpp");
    Rcpp::Function checks = base[".initialChecks"];
    // Call R checks and updata data and parameters
    return checks(Y0,Y1,P0,S0,params,diss,alpha,w);
    
  }catch (Rcpp::exception& e) {
    // Handle the Rcpp exception
    Rcpp::Rcerr << "Caught exception: " << e.what() << std::endl;
    return Rcpp::List::create();
    
  } catch (...){
    // Handle other types of exceptions
    Rcpp::Rcerr << "Caught unknown exception." << std::endl;
    return Rcpp::List::create();
  }
}

RCPP_EXPOSED_CLASS(ProbKMA);

RCPP_MODULE(ProbKMAModule) {
  Rcpp::class_<ProbKMA>("ProbKMA")
  .constructor<Rcpp::List,Rcpp::List,Rcpp::List,
               ProbKMA::matrix,ProbKMA::imatrix,std::string>()
  .method("probKMA_run",&ProbKMA::probKMA_run)
  .method("set_parameters", &ProbKMA::set_parameters);
}
