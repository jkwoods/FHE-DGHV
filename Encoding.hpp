//
//  Encoding.hpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//

#ifndef Encoding_hpp
#define Encoding_hpp

#include <stdio.h>
#include <vector>
#include "Pk.hpp"

class Encoding{
public:
    Encoding(Pk pk, std::vector<int> m); //for public creation
    ~Encoding();
    
    std::vector<int> decode();
    std::vector<int> decode_squashed();
    void recode();
    Encoding H_add(Encoding x);
    Encoding H_mult(Encoding x);
    
private:
    Encoding(Pk pk, mpz_class c); // class handling
    mpz_class e_val;
    Pk e_pk;

};

#endif /* Encoding_hpp */
