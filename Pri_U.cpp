//
//  Pri_U.cpp
//  FHE
//
//  Created by Woods, Jess on 10/17/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "Pri_U.hpp"
#include "Pk.hpp"
#include <cmath>
#include "utils.hpp"

Pri_U::Pri_U(Pk pk)
    : u_pk(pk), u_pri(makePri()), u_u(makeU()) {}

std::vector<int> Pri_U::getUList(){
    return u_u;
}

PseudoRandomInts Pri_U::makePri(){
    PseudoRandomInts pri = PseudoRandomInts(pow(2, u_pk.p_kap+1), u_pk.p_Theta);
    return pri;
}

std::vector<int> Pri_U::makeU(){
    std::vector<int> u = u_pri.getList(); //u draft
    
    for(int j = 0; j < u_pk.p_l; j++){
        std::vector<int> s1indices;
        int xpj = pow(2, u_pk.p_kap) / u_pk.p_p[j]; // i think its an int (??)
        
        std::vector<int> su;
        for(int i = 0; i < u_pk.p_Theta; i++){
            su.push_back(u_pk.p_s[j][i] * u[i]);
            if (u_pk.p_s[j][i] == 1){
                s1indices.push_back(i);
            }
        }
        
        int sumt = accumulate(su.begin(), su.end(), 0);
        sumt = mod(sumt, pow(2, u_pk.p_kap+1));
        
        while(sumt != xpj){
            //pick rand 1 in s
            int v = random_choice(s1indices);
            
            //change corresponding u
            su[v] = 0;
            int sumv = accumulate(su.begin(), su.end(), 0);
            int k1 = pow(2, u_pk.p_kap+1);
            int nu = k1 - sumv + xpj;
            
            while ((nu < 0) || (nu >= k1)){
                if (nu < 0){
                    nu = nu+k1;
                } else {
                    nu = nu-k1;
                }
            }
            
            u[v] = nu;
            
            //check not included - TODO if problems
            
        }
    }
    
    return u;
}

