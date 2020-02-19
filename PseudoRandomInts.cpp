//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "PseudoRandomInts.hpp"
#include "omp.h"

PseudoRandomInts::PseudoRandomInts(mpz_class x0, int len, long seed): r_x0(x0),  r_seed(seed), r_len(len), r_list(len)
{
    make_list();
}

PseudoRandomInts::PseudoRandomInts(mpz_class x0, int len): r_x0(x0),  r_seed(time(0)), r_len(len), r_list(len)
{
    make_list();
}

PseudoRandomInts::~PseudoRandomInts(){
    //TODO
}



void PseudoRandomInts::make_list(){
    //make class state
    gmp_randclass p_class_state (gmp_randinit_mt);
    p_class_state.seed(r_seed); // - now always the same (hopefully ... :)

//    #pragma omp parallel for
    for(int i = 0; i < r_len; i++){
        r_list[i] = p_class_state.get_z_range(r_x0); // 0 - (r_x0-1)
    }
}
