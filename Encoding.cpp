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
Encoding::Encoding(Pk pk, std::vector<int> m): e_pk(pk), e_val(pk.encode(m)) {} //public

Encoding::Encoding(Pk pk, mpz_class c): e_pk(pk), e_val(c) {} //private class handling

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

void Encoding::recode(int permutation_type){
    std::cout << "permuting" << "\n";
    e_val = e_pk.recode_and_permute(e_val);
}

Encoding Encoding::operator+(Encoding x){
    mpz_class c = e_pk.H_add(e_val, x.e_val);
    return Encoding(e_pk, c);
}

Encoding Encoding::operator*(Encoding x){
    mpz_class c = e_pk.H_mult(e_val, x.e_val);
    return Encoding(e_pk, c);
}

Encoding Encoding::neg(){
    Encoding one = Encoding(e_pk, {1,1,1,1,1,1,1,1,1,1});
    return *this+one;
}

Encoding Encoding::selector(std::vector<int> s, Encoding a, Encoding b){ // if s=1: a, else b

    Encoding sel = Encoding(a.e_pk, s);
    return (sel*a) + (sel.neg()*b);
}
