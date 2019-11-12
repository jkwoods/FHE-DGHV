//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "PseudoRandomInts.hpp"
#include "utils.hpp"

PseudoRandomInts::PseudoRandomInts(mpz_t x0, int len): r_len(len), r_seed(time(0)){
    
    make_r_x0(x0);
    makeList(); //state madestd::vector<mpz_t> r_list;

}

PseudoRandomInts::~PseudoRandomInts(){ //TODO - delete vectors??
    mpz_clear(r_x0);
    for(int i = 0; i < r_list.size(); i++) { mpz_clear(r_list[i]); }
}


long PseudoRandomInts::getSeed(){
    return r_seed;
}

void PseudoRandomInts::makeList(){
    
    gmp_randinit_mt(r_state);
    gmp_randseed_ui(r_state, r_seed);  //makeState - now always the same (hopefully ... :)
    
    for(int i = 0; i < r_len; i++){
        mpz_init(r_list[i]);
        mpz_urandomm(r_list[i], r_state, r_x0); //random_element(r_list[i], 0, r_x0);
    }
}

void PseudoRandomInts::make_r_x0(mpz_t x0){
    mpz_init_set(r_x0, x0);
}

