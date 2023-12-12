#include "ProbKMA.hpp"
#include "Factory.hpp"
#include <forward_list>
#include <limits>
#include <string_view>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// @TODO
// DARE CONSISTENZA A SET PARAMETERS PERCHE SE CAMBIO I PARAMETRI IN PROBKMA DEVONO
// ANCHE CAMBIARE IN MOTIFH1 ecc

class ProbKMA::_probKMAImp
{
public:
  
    _probKMAImp(const Rcpp::List& Y,const Rcpp::List& V,
                const Rcpp::List& parameters,
                const KMA::matrix& P0,const KMA::imatrix& S0,
                const std::string_view diss):
                _parameters(parameters),_P0(P0),_S0(S0)
               {  
                    // Initialize c++ Data structure
                    Initialize(Y,V,diss);
                    
                    // Create Dissimilarity factory
                    util::SharedFactory<Dissimilarity> dissfac;
                    dissfac.FactoryRegister<L2>("L2",_parameters._w);
                    dissfac.FactoryRegister<H1>("H1",_parameters._w,_parameters._alpha);
                    
                    //Create Motif factory
                    util::SharedFactory<MotifPure> motfac;
                    motfac.FactoryRegister<MotifL2>("L2");
                    motfac.FactoryRegister<MotifH1>("H1");
                    
                    //Create Performance factory
                    util::SharedFactory<PerformanceIndexAB> perfac;
                    perfac.FactoryRegister<PerformanceL2>("L2");
                    perfac.FactoryRegister<PerformanceH1>("H1");
                    
                    //check whether diss is valid
                    _motfac = motfac.instantiate(diss); // copy elision
                    _dissfac = dissfac.instantiate(diss); // copy elision
                    _perfac = perfac.instantiate(diss); // copy elision
                    if(not(_motfac and _dissfac and _perfac))
                      Rcpp::Rcerr<<"Invalid dissimilarity: Choose between L2,H1"<<std::endl;
                }

    ~_probKMAImp() = default;

    void Initialize(const Rcpp::List& Y,const Rcpp::List& V,std::string_view diss)
    {
      // Convert Rcpp Data Structure(List) into Armadillo data Structure(field)
      const Rcpp::List& Y0 = Y[0];
      const Rcpp::List& Y1 = Y[1];
      const Rcpp::List& V0 = V[0];
      const Rcpp::List& V1 = V[1];
      
      if (diss == "H1") {
        handleCaseH1(Y0, Y1, V0, V1);
      } else if (diss == "L2") {
        handleCaseL2(Y0, Y1, V0, V1);
      }
    }
      
      // Funzione di supporto per il caso "H1"
      void handleCaseH1(const Rcpp::List& Y0, const Rcpp::List& Y1,
                        const Rcpp::List& V0, const Rcpp::List& V1) {
        std::size_t Y_size = Y0.size();
        std::size_t V_size = V0.size();
        _Y.set_size(Y_size, 2);
        _V.set_size(V_size,2);
        
        for (int i = 0; i < Y_size; i++) {
          _Y(i, 0) = Rcpp::as<KMA::matrix>(Y0[i]);
          _Y(i, 1) = Rcpp::as<KMA::matrix>(Y1[i]);
        }
        
        for (int i = 0; i < V_size; i++) {
          _V(i, 0) = Rcpp::as<KMA::matrix>(V0[i]);
          _V(i, 1) = Rcpp::as<KMA::matrix>(V1[i]);
        }
      }
      
      // Funzione di supporto per il caso "L2"
      void handleCaseL2(const Rcpp::List& Y0, const Rcpp::List& Y1,
                        const Rcpp::List& V0, const Rcpp::List& V1) {
        if (!Rf_isNull(Y0[0])) {
          handleNonNullY(Y0, V0);
        } else {
          handleNonNullY(Y1, V1);
        }
      }
      
      // Funzione di supporto per il caso "L2" quando Y[0] non è nullo
      void handleNonNullY(const Rcpp::List& Y, const Rcpp::List& V) {
        std::size_t Y_size = Y.size();
        std::size_t V_size = V.size();
        _Y.set_size(Y_size, 1);
        _V.set_size(V_size, 1);
        
        for (int i = 0; i < Y_size; i++) {
          _Y(i, 0) = Rcpp::as<KMA::matrix>(Y[i]);
        }
        
        for (int i = 0; i < V_size; i++) {
          _V(i, 0) = Rcpp::as<KMA::matrix>(V[i]);
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
      const KMA::vector quantile4clean(_parameters._quantile4clean); // convert to vector to compute quantile
      KMA::matrix D;
      KMA::umatrix keep;
      std::size_t _n_rows_V = _V.n_rows;
      std::size_t _n_rows_Y = _Y.n_rows;
      std::vector<arma::urowvec> V_dom(_n_rows_V);
      KMA::Mfield V_new(_n_rows_V,_Y.n_cols); 

      while(iter < iter_max and BC_dist > _parameters._tol)
      {
        iter++;
        ///// clean motifs //////////////////////////
        const KMA::matrix P_old = _P0;
        const unsigned int& iter4clean = _parameters._iter4clean;
        const double& tol4clean = _parameters._tol4clean;
        if((iter>1)&&(!(iter%iter4clean))&&(BC_dist<tol4clean))
        {
          keep = D < arma::quantile(D,quantile4clean);
          const KMA::uvector empty_k = arma::find(arma::sum(keep,0));
          if(!empty_k.empty())
          {
            KMA::uvector index_min = arma::index_min(D,0);
            keep(empty_k(index_min),index_min).fill(1);
          }
          _P0.zeros();
          _P0.elem(arma::find(keep)).fill(1); // set one where values of keep are true
        }
        for(int i = 0;i < _n_rows_V;++i)
        {
          const arma::urowvec& V_dom_temp = util::findDomain<KMA::matrix>(_V(i,0));
          // capire se ha senso dichiararlo fuori
          const auto& V_new_variant = _motfac->compute_motif(V_dom_temp,_S0.col(i),
                                                            _P0.col(i),_Y,
                                                            _parameters._m);
          
          if(auto ptr_1 = std::get_if<MotifPure::indexField>(&V_new_variant))
          {
            const arma::sword& index = ptr_1->second;
            _S0.col(i) += index;
            V_new.row(i) = ptr_1->first; 
          }
          else
          {
            V_new.row(i) = *(std::get_if<KMA::Mfield>(&V_new_variant));
          }
          V_dom[i] = util::findDomain<KMA::matrix>(V_new(i,0));
        }
        
        if((iter>1)&&(!(iter%_parameters._iter4elong))&&(BC_dist<_parameters._tol4elong))
        {
          _motfac -> elongate_motifs(V_new,V_dom,_S0,_P0,
                                     _Y,D, _parameters,
                                     _perfac,_dissfac);
        }
        
////// find shift warping minimizing dissimilarities /////////////
        KMA::vector sd(2);
        arma::imat S_new(_n_rows_Y,_n_rows_V);
        arma::mat  D_new(_n_rows_Y,_n_rows_V);
        KMA::ivector c_k(_n_rows_V);
        const auto transform_function = [this](const KMA::matrix& V_new0) 
        {return std::floor(V_new0.n_rows * (1 - this->_parameters._max_gap));};
        const KMA::Mfield& V_new0 = V_new.col(0);
        std::transform(V_new0.begin(),V_new0.end(),c_k.begin(),transform_function);
        c_k.elem(c_k < _parameters._c) = _parameters._c;
        
#ifdef _OPENMP
#pragma omp parallel for collapse(2) firstprivate(sd)
#endif
        for (unsigned int i = 0; i < _n_rows_V; ++i)
          for (unsigned int j = 0; j < _n_rows_Y; ++j){ 
            sd = _dissfac->find_diss(_Y.row(j),V_new.row(i),_parameters._w,_parameters._alpha,
                                     c_k(i)); 
            S_new(j,i) = sd(0);
            D_new(j,i) = sd(1);
          }
        return Rcpp::List::create(Rcpp::Named("P0")=_P0,Rcpp::Named("S0")=_S0,
                                  Rcpp::Named("V_dom")=V_dom); 
      }
      return Rcpp::List::create();
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
    KMA::Mfield _Y;
    KMA::Mfield _V;
    
    //Motif and dissimilarity
    std::shared_ptr<MotifPure> _motfac;
    std::shared_ptr<Dissimilarity> _dissfac;
    std::shared_ptr<PerformanceIndexAB> _perfac;
    
    //Parameters
    Parameters _parameters;
    
    // Membership and shifting matrix
    KMA::matrix _P0;
    KMA::imatrix _S0;
    
};



///////// Implementation of funtions declared in the HEADER file ///////////////

ProbKMA::ProbKMA(const Rcpp::List& Y,const Rcpp::List& V,
                 const Rcpp::List& parameters,
                 const KMA::matrix& P0,const KMA::imatrix& S0,
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
               KMA::matrix,KMA::imatrix,std::string>()
  .method("probKMA_run",&ProbKMA::probKMA_run)
  .method("set_parameters", &ProbKMA::set_parameters);
}
