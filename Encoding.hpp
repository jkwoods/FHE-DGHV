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
    Encoding(mpz_t val, Pk pk);
    
    std::vector<int> decode();
    std::vector<int> decode_squashed();
    Encoding recode();
    Encoding H_add(Encoding x);
    Encoding H_mult(Encoding x);
private:
    mpz_t e_val;
    Pk e_pk;

};

#endif /* Encoding_hpp */
