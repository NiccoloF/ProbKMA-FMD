#include "Dissimilarity.hpp"

H1::H1(const KMA::vector& w,double alpha):SobolDiss(w),_alpha(alpha){};


double H1::computeDissimilarity(const KMA::Mfield& Y_i,
                                const KMA::Mfield& V_i) const
{
    return (1-_alpha) * this -> distance(Y_i(0,0),V_i(0,0)) + 
            _alpha * this -> distance(Y_i(0,1),V_i(0,1));
}


