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
    
    int getSeed();
    std::vector<mpz_t> getList();

private:
    mpz_t r_x0;
    int r_seed;
    std::vector<mpz_t> r_list;
    
    std::vector<mpz_t> makeList();
    
};

#endif /* PseudoRandomInts_hpp */
