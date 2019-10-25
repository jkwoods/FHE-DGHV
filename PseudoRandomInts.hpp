//
//  PseudoRandomInts.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef PseudoRandomInts_hpp
#define PseudoRandomInts_hpp

#include <stdio.h>
#include <vector>
#include <gmp.h>

class PseudoRandomInts{
public:
    PseudoRandomInts(mpz_t x0, int len);
    int r_len;
    std::vector<mpz_t> r_list;
    
    int getSeed();

private:
    mpz_t r_x0;
    int r_seed;
    
    void makeList();
    void make_r_x0(mpz_t x0);
    
};

#endif /* PseudoRandomInts_hpp */
