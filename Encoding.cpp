//
//  Encoding.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "Encoding.hpp"

//constructor
Encoding::Encoding(int val, Pk pk): e_val(val), e_pk(pk) {}

std::vector<int> Encoding::decode(){
    std::vector<int> m = e_pk.decode(e_val);
    return m;
}

std::vector<int> Encoding::decode_squashed(){
    std::vector<int> m = e_pk.decode_squashed(e_val);
    return m;
}

Encoding Encoding::recode(){
    int r = e_pk.recode(e_val);
    return Encoding(r,e_pk);
}

Encoding Encoding::H_add(Encoding x){
    int c = e_pk.H_add(e_val, x.e_val);
    return Encoding(c,e_pk);
}

Encoding Encoding::H_mult(Encoding x){
    int c = e_pk.H_mult(e_val, x.e_val);
    return Encoding(c,e_pk);
}
