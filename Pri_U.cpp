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

Pri_U::Pri_U(Pk pk)
    : u_pk(pk), u_pri(makePri()), u_u(pk.p_Theta){
        makeU();
    }

Pri_U::~Pri_U(){
    //TODO
}

PseudoRandomInts Pri_U::makePri(){
    mpz_class power = pow(2,u_pk.p_kap+1); //TODO pow
    PseudoRandomInts pri = PseudoRandomInts(power, u_pk.p_Theta);
    return pri;
}

void Pri_U::makeU(){
    u_u = u_pri.r_list; //u draft
    
    for(int j = 0; j < u_pk.p_l; j++){
        std::vector<int> s1indices;
        mpz_class xpj = pow(2, u_pk.p_kap) / u_pk.p_p[j]; // i think its an int (??) pow TODO
        
        std::vector<mpz_class> su(u_pk.p_Theta);
        for(int i = 0; i < u_pk.p_Theta; i++){
            su[i] = (u_pk.p_s[j][i] * u_u[i]);
            
            if (u_pk.p_s[j][i] == 1){
                s1indices.push_back(i);
            }
        }
        
        mpz_class sumt = accumulate(su.begin(), su.end(), 0); //replace all accumulates
        sumt = mod(sumt, pow(2, u_pk.p_kap+1));
        
        while(sumt != xpj){
            //pick rand 1 in s
            int v = random_choice(s1indices);
            
            //change corresponding u
            su[v] = 0;
            mpz_class sumv = accumulate(su.begin(), su.end(), 0);
            
            mpz_class k1 = pow(2, u_pk.p_kap+1);
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
                su[i] = u_u[i]*u_pk.p_s[j][i];
            }
            
            sumt = accumulate(su.begin(), su.end(), 0); //replace all accumulates
            sumt = mod(sumt, pow(2, u_pk.p_kap+1));
            
        }
    }
}

