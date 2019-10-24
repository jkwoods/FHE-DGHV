//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "Pk.hpp"
#include "Deltas.hpp"
#include "Pri_U.hpp"
#include <gmp.h>

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_n(n), p_q0(make_q0()), p_x(make_x()), p_xi(make_xi()), p_ii(make_ii()), p_B(Theta/theta), p_s(make_s()), p_vert_s(make_vert_s()), p_u(make_u()), p_y(make_y()), p_o(make_o()){
    
    make_p();
    make_pi();
    make_q0();
    make_x0();
    
}

int Pk::encode(std::vector<int> m){
    //m*xi
    std::vector<mpz_t> m_xi;
    for (int i = 0; i < p_l; i++){
        mpz_init(m_xi[i]);
        mpz_mul_si(m_xi[i], p_xi[i], m[i]);
        //m_xi.push_back(m[i]*p_xi[i]);
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

std::vector<int> Pk::decode(mpz_t c){ //TODO - mpz
    std::vector<int> m;
    for (int i = 0; i < p_l; i++){
        m.push_back(mod(modNear(c,p_p[i]),2));
    }
    return m;
}

std::vector<int> Pk::decode_squashed(mpz_t c){ //TODO gen
    std::vector<int> temp;
    return temp;
}

int Pk::recode(mpz_t r, mpz_t c){ //TODO gen
    return c;
}

void Pk::H_add(mpz_t added, mpz_t c1, mpz_t c2){ //TODO - mpz
    mpz_t pre_added;
    mpz_init(pre_added);
    mpz_add(pre_added, c1, c2);
    
    mpz_mod(added, pre_added, p_x0);
    
    mpz_clear(pre_added);
    //int c = mod(c1+c2,p_x0);
}

void Pk::H_mult(mpz_t multed, mpz_t c1, mpz_t c2){ //TODO - mpz
    mpz_t pre_mult;
    mpz_init(pre_mult);
    mpz_mul(pre_mult, c1, c2);
    
    mpz_mod(multed, pre_mult, p_x0);
    
    mpz_clear(pre_mult);
    
    //int c = mod(c1*c2,p_x0);
}



//private helper
void Pk::make_p(){
    //std::vector<int> p;
    for (int i = 0; i < p_l; i++){
        mpz_init(p_p[i]);
        mpz_ui_pow_ui(p_p[i], 2, p_eta-1);
        //p.push_back(random_prime(pow(2,p_eta-1), pow(2,p_eta)));
    }
}

void Pk::make_pi(){ //prod of all p[i]
    mpz_init_set_ui(p_pi,1); //int pi = 1;
    for (int i = 0; i < p_l; i++){
        mpz_mul(p_pi, p_pi, p_p[i]);
        //pi = pi*p_p[i];
    }
}

void Pk::make_q0(){
    mpz_init(p_q0);
    mpz_ui_pow_ui(p_q0, 2, p_gam);
    
    
    mpz_t comp;
    mpz_init(comp);
    mpz_ui_pow_ui(comp, 2, p_gam);
    mpz_fdiv_q(comp, comp, p_pi);

    //int q0 = pow(2,p_gam);
    
    mpz_t q0_prime1;
    mpz_init(q0_prime1);
    
    mpz_t q0_prime2;
    mpz_init(q0_prime2);
    
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, time(0)); //time as seed TODO - check this shit
    
    long lam2 = pow(p_lam,2);
    
    mpz_t nn;
    mpz_init(nn);
    mpz_ui_pow_ui(nn, 2, lam2);
    
    while (mpz_cmp(p_q0, comp)){     //while (q0 > (pow(2,p_gam)/p_pi)){
        
        //int q0_prime1 = random_prime(0, pow(2,pow(p_lam,2)));
        mpz_urandomm(q0_prime1, rstate, nn); // prime in 0, nn - 1
        while (mpz_probab_prime_p(q0_prime1, 50) == 0){ //isnt prime
            mpz_urandomm(q0_prime1, rstate, nn);
        }
        //int q0_prime2 = random_prime(0, pow(2,pow(p_lam,2)));
        mpz_urandomm(q0_prime2, rstate, nn);
        while (mpz_probab_prime_p(q0_prime2, 50) == 0){ //isnt prime
            mpz_urandomm(q0_prime2, rstate, nn);
        }
        
        mpz_mul(p_q0, q0_prime1, q0_prime2);
        //q0 = q0_prime1*q0_prime2;
    }
    mpz_clear(comp);
    mpz_clear(q0_prime1);
    mpz_clear(q0_prime2);
    mpz_clear(nn);

}

void Pk::make_x0(){
    mpz_init(p_x0);
    mpz_mul(p_x0, p_pi, p_q0);
}

std::vector<mpz_t> Pk::make_x(){
    Deltas x_D = Deltas(*this, p_tau, p_rhoi-1, 0);
    std::vector<mpz_t> x = x_D.getDeltaList();
    return x;
}

std::vector<mpz_t> Pk::make_xi(){
    Deltas xi_D = Deltas(*this, p_l, p_rho, 1);
    std::vector<mpz_t> xi = xi_D.getDeltaList();
    return xi;
}

std::vector<mpz_t> Pk::make_ii(){
    Deltas ii_D = Deltas(*this, p_l, p_rho, 2);
    std::vector<mpz_t> ii = ii_D.getDeltaList();
    return ii;
}

std::vector<std::vector<int>> Pk::make_s(){
    std::vector<std::vector<int>> s;
    for(int i = 0; i < p_l; i++){
        std::vector<int> sj;
        for (int j = 0; j < p_theta; j++){
            std::vector<int> fill;
            for (int b = 0; b < p_B; b++){
                fill.push_back(0);
            }
            if (j==0){
                fill[j] = 1;
                sj.insert(std::end(sj), std::begin(fill), std::end(fill));
            } else {
                sj.insert(std::end(sj), std::begin(fill), std::end(fill));
            }
        }
        s.push_back(sj);
    }
    
    for(int i = 0; i < p_theta; i++){ //TODO??
        std::vector<int> sri = random_sample(p_B, p_l);
        for(int j = 0; i < p_l; i++){
            int k = (p_B*i)+sri[j];
            s[j][k] = 1;
        }
    }
    return s;
}

std::vector<std::vector<int>> Pk::make_vert_s(){
    std::vector<std::vector<int>> vert_s;
    for(int i = 0; i < p_Theta; i++){
        std::vector<int> vsi;
        for(int j = 0; j < p_l; j++){
            vsi.push_back(p_s[j][i]);
        }
        vert_s.push_back(vsi);
    }
    return vert_s;
}

std::vector<mpz_t> Pk::make_u(){
    Pri_U priu = Pri_U(*this);
    std::vector<mpz_t> u = priu.getUList();
    return u;
}

std::vector<mpz_t> Pk::make_y(){
    std::vector<mpz_t> y;
    for (int i = 0; i < p_u.size(); i++){
        y.push_back(p_u[i] / pow(2, p_kap));
    }
    return y;
}

std::vector<mpz_t> Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    std::vector<mpz_t> o = o_D.getDeltaList();
    return o;
}

std::vector<int> Pk::random_sample(int range, int l){
    std::vector<int> sample;
    for(int i = 0; i < range; i++){
        sample.push_back(i);
    }
    std::random_shuffle(sample.begin(), sample.end());
    std::vector<int> cut_sample;
    for(int i = 0; i < l; i++){
        cut_sample.push_back(sample[i]);
    }
    return cut_sample;
}
