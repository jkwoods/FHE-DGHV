//
//  Deltas.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef Deltas_hpp
#define Deltas_hpp

#include <stdio.h>
#include "PseudoRandomInts.hpp"
#include "Pk.hpp"

class Deltas{
public:
    Deltas(Pk pk, int lenv, int rho, int cr);
    
    std::vector<int> getDeltaList();
    PseudoRandomInts r_pri;
    
private:
    
    std::vector<int> r_deltas;
    std::vector<int> r_Chi;
    Pk r_pk;
    int r_lenv;
    int r_rho;
    int r_cr;
    
    std::vector<int> makeDeltas();
    PseudoRandomInts makePri();
};

#endif /* Deltas_hpp */
