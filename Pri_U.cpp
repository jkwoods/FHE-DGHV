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
    : u_pk(pk), u_pri(makePri()){
        
        makeU();
    }

Pri_U::~Pri_U(){
    for(int i = 0; i < u_u.size(); i++) { mpz_clear(u_u[i]); }
}

PseudoRandomInts Pri_U::makePri(){
    mpz_t f_term;
    mpz_init(f_term);
    mpz_ui_pow_ui(f_term, 2, u_pk.p_kap+1);
    
    PseudoRandomInts pri = PseudoRandomInts(f_term, u_pk.p_Theta);
    
    mpz_clear(f_term);
    
    return pri;
}

void Pri_U::makeU(){
    //std::vector<int> u = u_pri.getList(); //u draft
    for (int i = 0; i < u_pri.r_len; i++){
        mpz_init_set(u_u[i], u_pri.r_list[i]);
    }
    
    mpz_t xpj;
    mpz_init(xpj);
    
    mpz_t two_kap;
    mpz_init(two_kap);
    mpz_ui_pow_ui(two_kap, 2, u_pk.p_kap);
    
    mpz_t two_kap_p1;
    mpz_init(two_kap_p1);
    mpz_ui_pow_ui(two_kap_p1, 2, u_pk.p_kap+1);
    
    std::vector<mpz_t> su;
    mpz_t sumsu;
    mpz_init(sumsu);
    
    mpz_t nu;
    mpz_init(nu);
    
    for(int j = 0; j < u_pk.p_l; j++){
        std::vector<int> s1indices;
        
        mpz_fdiv_q(xpj, two_kap, u_pk.p_p[j]);
        //int xpj = pow(2, u_pk.p_kap) / u_pk.p_p[j]; // i think its an int (??)
        
        //std::vector<int> su;
        for(int i = 0; i < u_pk.p_Theta; i++){
            mpz_init(su[i]);
            mpz_mul_ui(su[i], u_u[i], u_pk.p_s[j][i]); //su.push_back(u_pk.p_s[j][i] * u[i]);
            
            if (u_pk.p_s[j][i] == 1){
                s1indices.push_back(i);
            }
        }
        
        
        //int sumt = accumulate(su.begin(), su.end(), 0);
        for (int k = 0; k < su.size(); k++){
            mpz_add(sumsu, sumsu, su[k]);
        }
        
        mpz_mod(sumsu, sumsu, two_kap_p1); //sumt = mod(sumt, pow(2, u_pk.p_kap+1));
        
        while(mpz_cmp(sumsu, xpj)!=0){ //sumt != xpj
            //pick rand 1 in s
            int v = random_choice(s1indices);
            
            //change corresponding u
            mpz_set_ui(su[v], 0);
            
            //int sumv = accumulate(su.begin(), su.end(), 0);
            mpz_set_ui(sumsu, 0);
            for (int k = 0; k < su.size(); k++){
                mpz_add(sumsu, sumsu, su[k]);
            }
            
            //int k1 = pow(2, u_pk.p_kap+1);
            //int nu = k1 - sumv + xpj;
            mpz_sub(nu, two_kap_p1, sumsu);
            mpz_add(nu, nu, xpj);
            
            while ((mpz_cmp_ui(nu, 0)<0) || (mpz_cmp(nu, two_kap_p1)>=0)){ //((nu < 0) || (nu >= k1)){
                if (mpz_cmp_ui(nu, 0)<0){
                    mpz_add(nu, nu, two_kap_p1);//nu = nu+k1;
                } else {
                    mpz_sub(nu, nu, two_kap_p1); //nu = nu-k1;
                }
            }
            
            mpz_set(u_u[v], nu); //u[v] = nu;
            
            //su redo
            for(int i = 0; i < u_pk.p_Theta; i++){
                mpz_mul_ui(su[i], u_u[i], u_pk.p_s[j][i]);
            }
            
            mpz_set_ui(sumsu, 0);
            for (int k = 0; k < su.size(); k++){
                mpz_add(sumsu, sumsu, su[k]);
            }
            mpz_mod(sumsu, sumsu, two_kap_p1);
            
        }
    }
    
    //return u;
    

    for(int i = 0; i < u_pk.p_Theta; i++){ mpz_clear(su[i]); }
    mpz_clear(xpj);
    mpz_clear(two_kap);
    mpz_clear(two_kap_p1);
    mpz_clear(sumsu);
    mpz_clear(nu);
}

