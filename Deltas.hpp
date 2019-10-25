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

    PseudoRandomInts r_pri;
    std::vector<mpz_t> r_x;
    
private:
    
    std::vector<mpz_t> r_deltas;
    std::vector<mpz_t> r_Chi;
    
    Pk r_pk;
    int r_lenv;
    int r_rho;
    int r_cr;
    
    void makeDeltas();
    PseudoRandomInts makePri();
    void makeDeltaList();
};

#endif /* Deltas_hpp */
