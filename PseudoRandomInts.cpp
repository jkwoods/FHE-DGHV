//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "PseudoRandomInts.hpp"
#include "utils.hpp"

PseudoRandomInts::PseudoRandomInts(mpz_t x0, int len): r_len(len), r_seed(set_random_seed(0)){
    
    make_r_x0(x0);
    makeList();
    
}

PseudoRandomInts::~PseudoRandomInts(){ //TODO - delete vectors??
    mpz_clear(r_x0);
    for(int i = 0; i < r_list.size(); i++) { mpz_clear(r_list[i]); }
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

