//
//  Pri_U.cpp
//  FHE
//
//  Created by Woods, Jess on 10/17/19.
//  work done at Oak Ridge National Lab
//

#include "Pri_U.hpp"
#include "Pk.hpp"
#include <cmath>
#include "utils.hpp"

Pri_U::Pri_U(Pk& pk, int Theta)
: u_pk(pk), u_pri(makePri()), u_u(Theta)
{
    makeU();
}

Pri_U::~Pri_U(){
    //TODO
}

PseudoRandomInts Pri_U::makePri(){
    mpz_class pwr = power(2,u_pk.p_kap+1);
    PseudoRandomInts pri = PseudoRandomInts(pwr, u_pk.p_Theta);
    return pri;
}

void Pri_U::makeU(){
    for (int i = 0; i < u_pri.r_len; i++){
        u_u[i] = u_pri.r_list[i];
    } //u draft
    
    for(int j = 0; j < u_pk.p_l; j++){
        std::vector<int> s1indices;
        mpz_class xpj = power(2, u_pk.p_kap) / u_pk.p_p[j]; // i think its an int (??) pow TODO
        
        std::vector<mpz_class> su(u_pk.p_Theta);
        for(int i = 0; i < u_pk.p_Theta; i++){
            su[i] = (u_pk.p_s[j][i] * u_u[i]);
            
            if (u_pk.p_s[j][i] == 1){
                s1indices.push_back(i);
            }
        }
        
        mpz_class sumt = sum_array(su);
        sumt = sumt % power(2, u_pk.p_kap+1);
        
        while(sumt != xpj){
            //pick rand 1 in s
            int v = random_choice(s1indices);
            
            //change corresponding u
            su[v] = 0;
            mpz_class sumv = sum_array(su);
            
            mpz_class k1 = power(2, u_pk.p_kap+1);
            mpz_class nu = k1 - sumv + xpj;
            
            while ((nu < 0) || (nu >= k1)){
                if (nu < 0){
                    nu = nu+k1;
                } else {
                    nu = nu-k1;
                }
            }
            
            u_u[v] = nu;
            
            //su redo
            for(int i = 0; i < u_pk.p_Theta; i++){
                mpz_class temp = u_pk.p_s[j][i] * u_u[i];;
                su[i] = temp;
            }
            
            sumt = sum_array(su);
            sumt = sumt % power(2, u_pk.p_kap+1);
            
        }
    }
}

