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


void Encoding::recode(int shift=0){
    if (shift > 0){
        e_val = e_pk.recode_and_permute(e_val);
    } else {
        e_val = e_pk.recode(e_val);
    }
    
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

Encoding Encoding::matmul(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b, Pk pk){
    
    int row_a = a.size(); //first
    int col_a = a[0].size(); //second
    int row_b = b.size();
    int col_b = b[0].size();
    
    if (col_a != row_b){
        std::cout << "Error: Matrices cannot be multiplied.\n";
    }
    
    int mid_dim = col_a;
    
    int mults_per_row = mid_dim * col_b;
    int slots = row_a * mults_per_row;
    int plus = mid_dim;
    
    if (pk.p_l < slots){
        std::cout << "Error: This public key cannot encode these matrix sizes.\n";
    }
    std::vector<int> expanded_a(slots, 0);
    std::vector<int> expanded_b(slots, 0);
    
    //make expanded
    int e = 0;
    int ai = 0;
    int aj = 0;
    for (int t = 0; t < (row_a*col_a); t++){
        for (int r = 0; r < col_b; r++){
            
            expanded_a[e] = a[ai][aj];
            e++;
        }
        if (aj < col_a-1){
            aj++;
        } else {
            ai++;
            aj = 0;
        }
    }
    
    
    for (int r = 0; r < row_a; r++){
        e = r*col_b*row_b;
        int bi = 0;
        int bj = 0;
        std::cout << "e=" << e << "\n";
        for (int t = 0; t < (row_b*col_b); t++){
            std::cout << bi << "\n";
            std::cout << bj << "\n";
            expanded_b[e] = b[bi][bj];
            e++;
            if (bj < col_b-1){
                bj++;
            } else {
                bi++;
                bj = 0;
            }
        }
    }
    
    print_vec(expanded_a);
    print_vec(expanded_b);
    
    //make encoded
    Encoding encoded_a = Encoding(pk, expanded_a);
    Encoding encoded_b = Encoding(pk, expanded_b);
    
    //shift and add
    std::vector<Encoding> to_sum;
    Encoding result = encoded_a*encoded_b;
    for (int i = 1; i < plus; i++){
        to_sum.push_back(result);
        for (int j = 0; j < i; j++){ //shift correct number of times
            to_sum[i-1].recode(1);
        }
    }
    result.recode();
    
    for (int i = 0; i < plus-1; i++){
        result = result+to_sum[i];
    }
    
    return result;
}
