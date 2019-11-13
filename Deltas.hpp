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
public:
    Deltas(Pk pk, int lenv, int rho, int cr);
    ~Deltas();
    
    PseudoRandomInts r_pri;
    std::vector<mpz_class> r_x;
    
private:
    
    std::vector<mpz_class> r_deltas;
    std::vector<mpz_class> r_Chi;
    
    Pk r_pk;
    int r_lenv;
    int r_rho;
    int r_cr;
    gmp_randstate_t d_state; //for random crap - not psuedorandom PRI
    
    std::vector<mpz_class> makeDeltas();
    PseudoRandomInts makePri();
    std::vector<mpz_class> makeDeltaList();
    void makeState();
};

#endif /* Deltas_hpp */
