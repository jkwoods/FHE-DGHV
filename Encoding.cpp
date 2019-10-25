//
//  Encoding.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//

#include "Encoding.hpp"

//TODO - deconstructor - make sure you clear the e_val

//constructor
Encoding::Encoding(Pk pk, std::vector<int> m): e_pk(pk) {
    mpz_init(e_val);
    e_pk.encode(e_val, m);
}

//destructor
Encoding::~Encoding(){
    mpz_clear(e_val);
}

std::vector<int> Encoding::decode(){
    std::vector<int> m = e_pk.decode(e_val);
    return m;
}

std::vector<int> Encoding::decode_squashed(){
    std::vector<int> m = e_pk.decode_squashed(e_val);
    return m;
}

Encoding Encoding::recode(){
    mpz_t r;
    mpz_init(r);
    
    e_pk.recode(r, e_val);
    return Encoding(r,e_pk);
}

Encoding Encoding::H_add(Encoding x){
    mpz_t c;
    mpz_init(c);
    
    e_pk.H_add(c, e_val, x.e_val);
    return Encoding(c,e_pk);
}

Encoding Encoding::H_mult(Encoding x){
    mpz_t c;
    mpz_init(c);
    
    e_pk.H_mult(c, e_val, x.e_val);
    return Encoding(c,e_pk);
}
