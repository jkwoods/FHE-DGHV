//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "Pk.hpp"
#include <math.h>
#include <cmath>
#include <numeric>

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_n(n), p_p(make_p()), p_pi(make_pi()), p_q0(make_q0()), p_x0(p_pi*p_q0), p_x(make_x()), p_xi(make_xi()), p_ii(make_ii()), p_B(Theta/theta), p_s(make_s()), p_vert_s(make_vert_s()), p_u(make_u()), p_o(make_o()){}

int Pk::encode(std::vector<int> m){
    //m*xi
    std::vector<int> m_xi;
    for (int i = 0; i < p_l; i++){
        m_xi.push_back(m[i]*p_xi[i]);
    }
    
    //bi*ii
    std::vector<int> bi_ii;
    for (int i = 0; i < p_l; i++){
        int bi = (random_element(pow(-2,p_alphai),pow(2,p_alphai)));
        bi_ii.push_back(bi*p_ii[i]);
    }
    
    //b*x
    std::vector<int> b_x;
    for (int i = 0; i < p_tau; i++){
        int b = (random_element(pow(-2,p_alpha),pow(2,p_alpha)));
        b_x.push_back(b*p_x[i]);
    }
    
    int e_sum = accumulate(m_xi.begin(), m_xi.end(), 0) + accumulate(bi_ii.begin(), bi_ii.end(), 0) + accumulate(b_x.begin(), b_x.end(), 0);
    
    int c = modNear(e_sum,p_x0);
    
    return c;
}

std::vector<int> Pk::decode(int c){
    std::vector<int> m;
    for (int i = 0; i < p_l; i++){
        m.push_back(mod(modNear(c,p_p[i]),2));
    }
    return m;
}

//std::vector<int> Pk::decode_squashed(int c){}

//int Pk::recode(int c){}

//int Pk::H_add(int c1, int c2){}

//int Pk::H_mult(int c1, int c2){}

//private helper
//std::vector<int> Pk::make_p(){
    
//}

//int Pk::make_pi(){
    
//}

//int Pk::make_q0(){
    
//}

//std::vector<int> Pk::make_x(){
    
//}

//std::vector<int> Pk::make_xi(){
    
//}

//std::vector<int> Pk::make_ii(){
    
//}

//std::vector<std::vector<int>> Pk::make_s(){
    
//}

//std::vector<std::vector<int>> Pk::make_vert_s(){
    
//}

//std::vector<int> Pk::make_u(){
    
//}

//std::vector<int> Pk::make_o(){
    
//}


