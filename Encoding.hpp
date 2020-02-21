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

    Encoding(Pk pk, mpz_class c); // class handling

    ~Encoding();
    
    std::vector<int> decode();
    std::vector<int> decode_squashed();
    void recode(int shift);
    Encoding operator+(Encoding x);
    Encoding operator*(Encoding x);
    static Encoding selector(std::vector<int> s, Encoding a, Encoding b);
    Encoding neg();
    mpz_class e_val;
    
    
private:
    Pk e_pk;
    
};

#endif /* Encoding_hpp */
