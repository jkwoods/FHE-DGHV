//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "PseudoRandomInts.hpp"


PseudoRandomInts::PseudoRandomInts(mpz_class x0, int len): r_x0(x0), r_len(len), r_seed(time(0)), r_list(make_list()){
}

PseudoRandomInts::~PseudoRandomInts(){
    //TODO
}

long PseudoRandomInts::getSeed(){
    return r_seed;
}

std::vector<mpz_class> PseudoRandomInts::make_list(){
    
    gmp_randinit_mt(r_state);
    gmp_randseed_ui(r_state, r_seed);  //makeState - now always the same (hopefully ... :)
    
    std::vector<mpz_class> list(r_len);
    for(int i = 0; i < r_len; i++){
        list[i] = random_element(0, r_x0);
    }
    
    return list;
}
