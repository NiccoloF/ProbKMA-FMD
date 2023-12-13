#include "ProbKMA.hpp"
#include "Factory.hpp"
#include <forward_list>
#include <limits>
#include <string_view>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp20)]]

// @TODO: dare consistenza a set parameters, perchè se cambio i parametri in probKma devo anche cambiarli in MotifH1 ecc

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
      const unsigned int& iter_max = _parameters._iter_max;
      // Riccardo comment: why std::forward_list ? 
      KMA::vector J_iter(iter_max,arma::fill::zeros);
      KMA::vector BC_dist_iter(iter_max,arma::fill::zeros);
      auto BC_dist = std::numeric_limits<double>::infinity();
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
        KMA::imatrix S_new(_n_rows_Y,_n_rows_V);
        KMA::matrix  D_new(_n_rows_Y,_n_rows_V);
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

        // compute memberships (too much code in the run?!)
        // @TODO: change types to KMA:: ...
        KMA::matrix P_new(_n_rows_Y,_n_rows_V,arma::fill::zeros);
        KMA::umatrix D0 = (D_new == 0);
        const KMA::uvector & mult_assign = arma::find(arma::sum(D0,1) > 1);
        for (arma::sword i : mult_assign) {
          // @TODO: complete this warning message as the message of the prof
          Rcpp::warning("Curve has dissimilarity 0 from two different motifs. Using only one of them..."); 
          const KMA::uvector & indexes = arma::find(D0.row(i) == 1);
          D0.row(i).zeros();
          D0(i,indexes(arma::randi(arma::distr_param(0, indexes.n_elem - 1)))) = 1;
        }
        const KMA::uvector & D0_index = arma::find(arma::sum(D0,1) == 1);
        for(arma::sword i : D0_index) {
          const KMA::uvector & col = arma::find(D0.row(i)==1);
          P_new(i,col(0)) = 1;
        }
        const KMA::uvector & not_D0_index = arma::find(arma::sum(D0,1) !=1);
        const KMA::matrix & Dm = arma::pow(D_new.rows(not_D0_index),1/(_parameters._m-1));
        P_new.rows(not_D0_index) = 1 / (Dm % arma::repmat(arma::sum(1/Dm,1),1,Dm.n_cols));
        const KMA::uvector & deg_indexes = arma::find(arma::sum(P_new,0)==0);
        for (arma::sword k : deg_indexes) {
          // @TODO: complete this warning message as the message of the prof
          Rcpp::warning("Motif is degenerate (zero membership). Selecting a new center..."); 
          P_new(index_min(D_new.col(k)),k) = 1;
        }

        // evaluate objective functions
        // Riccardo: new part, just written 
        KMA::matrix temp_DP = D_new % (arma::pow(P_new,_parameters._m));
        temp_DP.replace(arma::datum::nan,0);
        J_iter(iter) = arma::accu(temp_DP);
      
        // compute Bhattacharyya distance between P_old and P_new
        // @TODO: ask to professor why RowSums in her code and not ColSums
        const arma::rowvec & BC_dist_k = -arma::log(arma::sum(arma::sqrt(P_old % P_new),0));
        std::string_view criterion = _parameters._stopCriterion;
        if (criterion == "max")
          BC_dist = arma::max(BC_dist_k);
        else if (criterion == "mean")
          BC_dist = arma::mean(BC_dist_k);
        else if (criterion == "quantile")
          BC_dist = arma::conv_to<arma::vec>::from
                    (arma::quantile(BC_dist_k,arma::vec(_parameters._prob)))(0);
  
        BC_dist_iter(iter) = BC_dist;

        // update 
        _V = V_new;
        _P0 = P_new; 
        _S0 = S_new;
        D = D_new; 
      }

      /////  prepare output //////////////////////////////////
      KMA::matrix  P_clean(_n_rows_V,_n_rows_Y,arma::fill::zeros);
      KMA::imatrix S_clean(_S0);
      KMA::matrix  D_clean(_n_rows_Y,_n_rows_V);
      // Riccardo comment: la prof non fa nessun controllo in questa parte su quali delle due strutture viene restituita, perchè?
      // non viene fatto controllo se pair_motif contiene effettivamente solo il campo arma::field<arma::mat>
      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto & pair_motif_shift = _motfac->compute_motif(V_dom[k], _S0.col(k),
                                                              _P0.col(k), _Y,
                                                              _parameters._m);
        _V.row(k) = *(std::get_if<arma::field<arma::mat>>(&pair_motif_shift));
      } 
      KMA::umatrix keep = D < arma::quantile(D,quantile4clean);
      KMA::uvector empty_k = arma::find(arma::sum(keep,0) == 0);
      for (arma::sword k: empty_k)
        keep(arma::index_min(D.col(k)),k) = 1;
      P_clean(arma::find(keep)) = 1;
      KMA::Mfield V_clean(_V.n_cols,_n_rows_V); // come impostare la size giusta di V_clean? usando _V.n_cols?
      std::map<arma::sword,arma::sword> shift_s;
      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto & new_motif =  _motfac->compute_motif(V_dom[k], _S0.col(k),
                                                         P_clean.col(k), _Y,
                                                         _parameters._m);
        if (auto ptr = std::get_if<arma::field<arma::mat>>(&new_motif)){
          V_clean.col(k) = *ptr;
        } else {
          const auto & pair_motif_shift = *(std::get_if<std::pair<arma::field<arma::mat>,arma::sword>>(&new_motif));
          V_clean.col(k) = pair_motif_shift.first;
          shift_s.insert(std::make_pair(k, pair_motif_shift.second));
        }
      }
      for(auto it = shift_s.begin();it != shift_s.cend(); ++it){
        S_clean.col(it->first) += it->second;
      }

      std::vector<KMA::uvector> V_dom_new(K); // questi vector di urowvec -> field<urowvec> per consistenza?
      for(arma::uword k=0; k < _n_rows_V ; ++k){
        V_dom_new[k] = util::findDomain<KMA::matrix>(_V(k,0));
      }
      // compute dissimilarities from cleaned motifs
      // Riccardo comment: da rivedere con molta cura: COMPUTATIONAL COST + TEMPLATE VERSION per use0,use1 + dichiarare fuori?
      const arma::uword d = Y(0,0).n_cols;
      arma::uword index_row;
      arma::uword index_size;
      KMA::matrix y0;
      KMA::matrix y1;
      KMA::Mfield y(1,_Y.n_cols); // in this way should be 1 or 2 according to use0, use1
      for(arma::uword k=0; k < _n_rows_V; ++k){
        const auto & s_k = S.col(k);
        const auto & v_dom = V_dom[k];
        const int v_len = v_dom.size();
        KMA::Mfield v_clean = V_clean.col(k);
        const KMA::uvector & indeces_dom = arma::find(v_dom==0);
        for (arma::uword i=0; i < _n_rows_Y; ++i){
          const int s = s_k(i);
          KMA::ivector index = std::max(1,s) - 1 + arma::regspace<arma::ivec>(1,v_len - std::max(1,s));
          index_size = index.size();
          v_clean(0,k).shed_rows(indeces_dom);
          const arma::uword y_len = Y(i,0).n_rows;
          y0.set_size(index_size + std::max(1,s),d);
          y0.fill(arma::datum::nan);
          for(unsigned int j = 0; j < index_size; ++j) {
            if (index[j]  <= y_len){
              index_row = std::max(0, 1-s) + j;
              y0.row(index_row) =  Y(i,0).row(index[j] - 1);
            }
          }
          y0.shed_rows(indeces_dom);
          y(0,0) = y0;
          if (use1){
            v_clean(1,k).shed_rows(indeces_dom);
            const arma::uword y_len = Y(i,1).n_rows;
            y1.set_size(index_size + std::max(1,s),d);
            y1.fill(arma::datum::nan);
            for(unsigned int j = 0; j < index_size; ++j) {
              if (index[j]  <= y_len){
                index_row = std::max(0, 1-s) + j;
                y1.row(index_row) =  Y(i,1).row(index[j] - 1);
              }
            }
            y1.shed_rows(indeces_dom);
            y(0,1) = y1;
          }
        D_clean(i,k) = _dissfac->computeDissimilarity(y,V_clean.col(i)); 
        }
      }

      return toR(V_clean,P_clean,S_clean,D,D_clean,J_iter,BC_dist_iter,iter);
    }

  
    void set_parameters(const Rcpp::List& newParameters)
    {
      _parameters = newParameters;
    }
    
    
    // return a R structure, V, _P0, _S0 già presenti <- convertire tutte le cose da convertire V, V_clean
    Rcpp::List toR(const KMA::Mfield & V_clean,
                   const KMA::matrix & P_clean,
                   const KMA::imatrix & S_clean,
                   const KMA::matrix & D,
                   const KMA::matrix & D_clean,
                   const KMA::vector & J_iter,
                   const KMA::vector & BC_dist_iter
                   std::size_t iter) const  
    {   
        // conv to Rcpp::List of V0,V1,V0_clean,V1_clean
        Rcpp::List V0(_V.n_rows);
        Rcpp::List V1(_V.n_rows);
        Rcpp::List V0_clean(_V.n_rows);
        Rcpp::List V1_clean(_V.n_rows);
        // Riccardo comment: devo capire in qualche modo come distinguere caso solo derivate o solo curve
        // Riccardo comment: una soluzione potrebbe essere salvarsi std::string_view diss che viene passata al constructor
        // Riccardo comment: per ora fillo solo V0 in entrambi casi ma @TODO: distinguere quando fare fill di V0 oppure fill di V1
        arma::uword use1 = _V.n_cols == 2? 1 : 0;
        for (arma::uword k = 0; k < _V.n_rows; ++k){
          V0[k] = _V(k,0);
          V0_clean = V_clean(0,k);
          if (use1) {
            V1[k] = _V(k,1);
            V1_clean = V_clean(1,k);
          }
        } 
        if (!(_parameters._return_options)){
          return Rcpp::List::create(Rcpp::Named("V0") = V0,
                                    Rcpp::Named("V1") = V1,
                                    Rcpp::Named("V0_clean") = V0_clean,
                                    Rcpp::Named("V1_clean") = V1_clean,
                                    Rcpp::Named("P0") = _P0,
                                    Rcpp::Named("P_clean") = P_clean,
                                    Rcpp::Named("S0") = _S0,
                                    Rcpp::Named("S_clean") = S_clean,
                                    Rcpp::Named("D") = D,
                                    Rcpp::Named("D_clean") = D_clean,
                                    Rcpp::Named("iter") = iter,
                                    Rcpp::Named("J_iter") = J_iter,
                                    Rcpp::Named("BC_dist_iter") = BC_dist_iter);
        } else {
          return Rcpp::List::create(Rcpp::Named("V0") = V0,
                                    Rcpp::Named("V1") = V1,
                                    Rcpp::Named("V0_clean") = V0_clean,
                                    Rcpp::Named("V1_clean") = V1_clean,
                                    Rcpp::Named("P0") = _P0,
                                    Rcpp::Named("P_clean") = P_clean,
                                    Rcpp::Named("S0") = _S0,
                                    Rcpp::Named("S_clean") = S_clean,
                                    Rcpp::Named("D") = D,
                                    Rcpp::Named("D_clean") = D_clean,
                                    Rcpp::Named("iter") = iter,
                                    Rcpp::Named("J_iter") = J_iter,
                                    Rcpp::Named("BC_dist_iter") = BC_dist_iter,
                                    Rcpp::Named("standardize") = _parameters._standardize,
                                    Rcpp::Named("K") = _parameters._K,
                                    Rcpp::Named("c") = _parameters._c,
                                    Rcpp::Named("c_max") = _parameters._c_max,
                                    Rcpp::Named("iter_max") = _parameters._iter_max,
                                    Rcpp::Named("quantile") = _parameters._quantile,
                                    Rcpp::Named("tol") = _parameters._tol,
                                    Rcpp::Named("stop_criterion") = _parameters._stopCriterion,
                                    Rcpp::Named("m") = _parameters._m,
                                    Rcpp::Named("iter4elong") = _parameters._iter4elong,
                                    Rcpp::Named("tol4elong") = _parameters._tol4elong,
                                    Rcpp::Named("max_elong") = _parameters._max_elong,
                                    Rcpp::Named("trials_elong") = _parameters._trials_elong,
                                    Rcpp::Named("deltaJk_elong") = _parameters.__deltaJK_elong,
                                    Rcpp::Named("max_gap") = _parameters._max_gap,
                                    Rcpp::Named("iter4clean") = _parameters._iter4clean,
                                    Rcpp::Named("tol4clean") = _parameters._tol4clean);
        }
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
