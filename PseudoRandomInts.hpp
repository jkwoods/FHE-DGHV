//
//  PseudoRandomInts.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#ifndef PseudoRandomInts_hpp
#define PseudoRandomInts_hpp

#include <stdio.h>
#include <vector>
#include <gmp.h>

class PseudoRandomInts{
public:
    PseudoRandomInts(mpz_t x0, int len);
    ~PseudoRandomInts();
    
    int r_len;
    int r_zero_temp;
    std::vector<mpz_t> r_list;
    
    long getSeed();

private:
    mpz_t r_x0;
    long r_seed;
    gmp_randstate_t r_state;
    
    void makeList();
    void make_r_x0(mpz_t x0);
    
};

#endif /* PseudoRandomInts_hpp */
