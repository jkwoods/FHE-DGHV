//
//  mpi_utils.hpp
//  FHE
//
//  Created by Woods, Jess on 1/16/20.
//  Copyright © 2020 Woods, Jess. All rights reserved.
//

#ifndef MultEncodings_hpp
#define MultEncodings_hpp

#include <stdio.h>
#include "Encoding.hpp"

class MultEncodings{
public:
    MultEncodings(Pk pk, std::vector<std::vector<int>> m_array);
    ~MultEncodings();
    
    std::vector<mpz_class> m_val_array;
    
    std::vector<std::vector<int>> decode();
    void recode();
    MultEncodings operator+(MultEncodings x);
    MultEncodings operator*(MultEncodings x);
    void shift(int offset);
    Encoding get_Encoding(int index);
    
private:
    MultEncodings(Pk pk, std::vector<mpz_class> c_array);
    
    Pk m_pk;
    
    void make_encodings(std::vector<std::vector<int>> m_array);

};

#endif /* MultEncodings_hpp */
