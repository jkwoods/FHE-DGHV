//
//  mpi_utils.cpp
//  FHE
//
//  Created by Woods, Jess on 1/16/20.
//  Copyright © 2020 Woods, Jess. All rights reserved.
//

#include "MultEncodings.hpp"

MultEncodings::MultEncodings(Pk pk, std::vector<std::vector<int>> m_array)
: m_pk(pk), m_val_array(m_array.size())
{
    make_encodings(m_array);
}

MultEncodings::MultEncodings(Pk pk, std::vector<mpz_class> c_array) //private
: m_pk(pk), m_val_array(c_array)
{}


MultEncodings::~MultEncodings(){
    //TODO -vectors?
}

void MultEncodings::make_encodings(std::vector<std::vector<int>> m_array){
    //TODO: MPI here
    
    for(int i = 0; i < m_array.size(); i++){
        m_val_array[i] = m_pk.encode(m_array[i]);
    }
}

std::vector<std::vector<int>> MultEncodings::decode(){
    //TODO: MPI here
    
    std::vector<std::vector<int>> m_result(m_val_array.size(),std::vector<int>(m_pk.p_l));
    
    for(int i = 0; i < m_val_array.size(); i++){
        m_result[i] = m_pk.decode(m_val_array[i]);
    }
    return m_result;

}

void MultEncodings::recode(){
    //TODO: MPI here
    
    for(int i = 0; i < m_val_array.size(); i++){
        m_val_array[i] = m_pk.recode(m_val_array[i]);
    }
    
}

MultEncodings MultEncodings::operator+(MultEncodings x){
    if (m_val_array.size() >= x.m_val_array.size()){
        std::vector<mpz_class> sum(x.m_val_array.size());
        
        for (int i = 0; i < x.m_val_array.size(); i++){
            sum[i] = m_pk.H_add(m_val_array[i], x.m_val_array[i]);
        }
        return MultEncodings(m_pk, sum);
    } else {
        std::vector<mpz_class> sum(m_val_array.size());
        
        for (int i = 0; i < m_val_array.size(); i++){
            sum[i] = m_pk.H_add(m_val_array[i], x.m_val_array[i]);
        }
        return MultEncodings(m_pk, sum);
    }
    
}

MultEncodings MultEncodings::operator*(MultEncodings x){
    if (m_val_array.size() >= x.m_val_array.size()){
        std::vector<mpz_class> prod(x.m_val_array.size());
        
        for (int i = 0; i < x.m_val_array.size(); i++){
            prod[i] = m_pk.H_mult(m_val_array[i], x.m_val_array[i]);
        }
        return MultEncodings(m_pk, prod);
    } else {
        std::vector<mpz_class> prod(m_val_array.size());
        
        for (int i = 0; i < m_val_array.size(); i++){
            prod[i] = m_pk.H_mult(m_val_array[i], x.m_val_array[i]);
        }
        return MultEncodings(m_pk, prod);
    }
    
}

void MultEncodings::shift(int offset){ //only takes positive TODO
    std::vector<mpz_class> new_val_array(m_val_array.size());
    for (int i = offset; i < m_val_array.size(); i++){
        new_val_array[i-offset] = m_val_array[i];

    }
    for (int i = 0; i < offset; i++){
        new_val_array[m_val_array.size()-offset+i] = m_val_array[i];
    }

    m_val_array = new_val_array;
}

Encoding MultEncodings::get_Encoding(int index){
    Encoding c = Encoding(m_pk, m_val_array[index]);
    return c;
}
