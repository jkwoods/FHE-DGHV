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
#include <gmpxx.h>
#include "utils.hpp"

class PseudoRandomInts{
private:
    mpz_class r_x0;
    long r_seed;
    void make_list();
    
public:
    PseudoRandomInts(mpz_class x0, int len, long seed);
    PseudoRandomInts(mpz_class x0, int len);
    ~PseudoRandomInts();
    
    int r_len;
    std::vector<mpz_class> r_list;
    


    
};

#endif /* PseudoRandomInts_hpp */
