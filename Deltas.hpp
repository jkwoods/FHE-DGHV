//
//  Deltas.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#ifndef Deltas_hpp
#define Deltas_hpp

#include <stdio.h>
#include "PseudoRandomInts.hpp"
#include "Pk.hpp"
#include "utils.hpp"
#include <cmath>

class Deltas{
    
private:
    Pk r_pk;
    int r_rho;
    int r_cr;
    
    std::vector<mpz_class> r_deltas;
    std::vector<mpz_class> r_Chi;
    
    void makeDeltas();
    void makeDeltaList();

public:
    Deltas(Pk& pk, int lenv, int rho, int cr);
    ~Deltas();
    
    
    int r_lenv;
    PseudoRandomInts r_pri;
    std::vector<mpz_class> r_x;
    
};

#endif /* Deltas_hpp */
