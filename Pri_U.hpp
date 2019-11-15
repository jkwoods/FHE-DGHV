//
//  Pri_U.hpp
//  FHE
//
//  Created by Woods, Jess on 10/17/19.
//  work done at Oak Ridge National Lab
//

#ifndef Pri_U_hpp
#define Pri_U_hpp

#include <stdio.h>
#include <vector>
#include "Pk.hpp"
#include "PseudoRandomInts.hpp"

class Pri_U{

    
private:
    
    Pk u_pk;
    
    void makeU();
    PseudoRandomInts makePri();
    
public:
    Pri_U(Pk& pk, int Theta);
    ~Pri_U();
    
    PseudoRandomInts u_pri;
    std::vector<mpz_class> u_u;
};


#endif /* Pri_U_hpp */
