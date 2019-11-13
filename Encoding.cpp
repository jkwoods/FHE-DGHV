//
//  Encoding.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//

#include "Encoding.hpp"

//TODO - deconstructor - make sure you clear the e_val

//constructors
Encoding::Encoding(Pk pk, std::vector<int> m): e_pk(pk), e_val(pk.encode(m)) {}

Encoding::Encoding(Pk pk, mpz_class c): e_pk(pk), e_val(c) {}

//destructor
Encoding::~Encoding(){
    //TODO
}

std::vector<int> Encoding::decode(){
    std::vector<int> m = e_pk.decode(e_val);
    return m;
}

std::vector<int> Encoding::decode_squashed(){
    std::vector<int> m = e_pk.decode_squashed(e_val);
    return m;
}

void Encoding::recode(){
    e_val = e_pk.recode(e_val);
}

Encoding Encoding::H_add(Encoding x){
    mpz_class c = e_pk.H_add(e_val, x.e_val);
    return Encoding(e_pk, c);
}

Encoding Encoding::H_mult(Encoding x){
    mpz_class c = e_pk.H_mult(e_val, x.e_val);
    return Encoding(e_pk, c);
}
