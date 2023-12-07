#include "Motif.hpp"

std::variant<MotifPure::indexField,KMA::Mfield>
Motif_L2::compute_motif(const arma::urowvec& v_dom,
                        const KMA::ivector& s_k,
                        const KMA::vector& p_k,
                        const KMA::Mfield& Y,
                        double m) const
  {
    return compute_motif_helper<false>(v_dom,s_k,p_k,Y,m);
  }
    
    
    