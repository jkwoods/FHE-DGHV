//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "PseudoRandomInts.hpp"
#include "utils.hpp"

PseudoRandomInts::PseudoRandomInts(mpz_t x0, int len): r_len(len), r_seed(set_random_seed(0)){
    
    make_r_x0(x0);
    makeList();
    
}

int PseudoRandomInts::getSeed(){
    return r_seed;
}

void PseudoRandomInts::makeList(){
    
    set_random_seed(r_seed);
    for(int i = 0; i < r_len; i++){
        mpz_init(r_list[i]);
        
        random_element(r_list[i], 0, r_x0); //not power of 2
    }
}

void PseudoRandomInts::make_r_x0(mpz_t x0){
    mpz_init_set(r_x0, x0);
}

