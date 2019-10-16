//
//  PseudoRandomInts.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "PseudoRandomInts.hpp"
#include "utils.hpp"

PseudoRandomInts::PseudoRandomInts(int x0, int len): r_x0(x0), r_len(len), r_seed(set_random_seed(0)), r_list(makeList()){}

int PseudoRandomInts::getSeed(){
    return r_seed;
}

std::vector<int> PseudoRandomInts::getList(){
    return r_list;
}

std::vector<int> PseudoRandomInts::makeList(){
    set_random_seed(r_seed);
    std::vector<int> prl;
    for(int i = 0; i < r_len; i++){
        int a = random_element(0, r_x0);
        prl.push_back(a);
    }
    return prl;
}

