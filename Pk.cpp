//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

// Important note to those trying to read code: gmp library (all the mpz_... stuff) is call by reference

#include "Pk.hpp"
#include "Deltas.hpp"
#include "Pri_U.hpp"
#include "Encoding.hpp"

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_n(n), p_B(Theta/theta), p_s(make_s()), p_vert_s(make_vert_s()) {
    
    make_p();
    make_pi();
    make_q0();
    make_x0();
    make_x();
    make_xi();
    make_ii();
    make_u();
    make_y();
    make_o();
    
}

void Pk::encode(mpz_t c, std::vector<int> m){
    //m*xi
    mpz_t m_xi;
    mpz_init(m_xi);
    for (int i = 0; i < p_l; i++){
        mpz_t m_xi_temp;
        mpz_init(m_xi_temp);
        mpz_mul_ui(m_xi_temp, p_xi[i], m[i]);
        
        //m_xi.push_back(m[i]*p_xi[i]);
        
        mpz_add(m_xi, m_xi, m_xi_temp);
        mpz_clear(m_xi_temp);
        
    }
    
    //bi*ii
    mpz_t bi_ii;
    mpz_init(bi_ii);
    for (int i = 0; i < p_l; i++){
        mpz_t bi_ii_temp;
        mpz_init(bi_ii_temp);
        
        mpz_t bi;
        mpz_init(bi);
        random_element(bi, p_alphai, p_alphai); //needs to raise these to -2 and 2
        
        mpz_mul(bi_ii_temp, p_ii[i], bi);
        
        mpz_add(bi_ii, bi_ii, bi_ii_temp);
        
        mpz_clear(bi);
        mpz_clear(bi_ii_temp);
        
        //int bi = (random_element(pow(-2,p_alphai),pow(2,p_alphai)));
        //bi_ii.push_back(bi*p_ii[i]);
    }
    
    //b*x
    mpz_t b_x;
    mpz_init(b_x);
    for (int i = 0; i < p_tau; i++){
        mpz_t b_x_temp;
        mpz_init(b_x_temp);
        
        mpz_t b;
        mpz_init(b);
        random_element(b, p_alpha, p_alpha); //needs to raise these to -2 and 2
        
        mpz_mul(b_x, p_x[i], b);
        mpz_add(b_x, b_x, b_x_temp);
        
        mpz_clear(b);
        mpz_clear(b_x_temp);
        
        //int b = (random_element(pow(-2,p_alpha),pow(2,p_alpha)));
        //b_x.push_back(b*p_x[i]);
    }
    
    //summation
    mpz_t bigsum;
    mpz_init(bigsum);
    mpz_add(bigsum, m_xi, bi_ii);
    mpz_add(bigsum, bigsum, b_x);
    
    modNear(c, bigsum, p_x0); //TODO

}

std::vector<int> Pk::decode(mpz_t c){
    std::vector<int> m;
    
    mpz_t b;
    mpz_init(b);
    
    mpz_t one;
    mpz_init_set_ui(one, 1);
    for (int i = 0; i < p_l; i++){
        modNear(b,c,p_p[i]);
        mpz_and(b, b, one);
        
        int last_bit = int (mpz_get_d(b));
        m.push_back(last_bit);
        
        //m.push_back(mod(modNear(c,p_p[i]),2));
    }
    
    mpz_clear(b);
    mpz_clear(one);
    
    return m;
}

std::vector<int> Pk::decode_squashed(mpz_t c){ //TODO gen
    std::vector<int> temp;
    return temp;
}

void Pk::recode(mpz_t r, mpz_t c){ //TODO gen
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

void Pk::make_x(){ //TODO DELTAS - TODO initialize list
    Deltas x_D = Deltas(*this, p_tau, p_rhoi-1, 0);
    for (int i = 0; i < p_x.size(); i++){
        mpz_init_set(p_x[i], x_D.r_x[i]);
    }
    //std::vector<mpz_t> x = x_D.getDeltaList();
    //return x;
}

void Pk::make_xi(){
    Deltas xi_D = Deltas(*this, p_l, p_rho, 1);
    for (int i = 0; i < p_xi.size(); i++){
        mpz_init_set(p_xi[i], xi_D.r_x[i]);
    }
    //std::vector<mpz_t> xi = xi_D.getDeltaList();
    //return xi;
}

void Pk::make_ii(){
    Deltas ii_D = Deltas(*this, p_l, p_rho, 2);
    for (int i = 0; i < p_ii.size(); i++){
        mpz_init_set(p_ii[i], ii_D.r_x[i]);
    }
    //std::vector<mpz_t> ii = ii_D.getDeltaList();
    //return ii;
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

void Pk::make_u(){
    Pri_U priu = Pri_U(*this);
    for (int i = 0; i < priu.u_u.size(); i++){
         mpz_init_set(p_u[i], priu.u_u[i]);
    }
    //std::vector<mpz_t> u = priu.getUList();
    //return u;
}

void Pk::make_y(){
    //std::vector<mpz_t> y;
    mpz_t two_kap;
    mpz_init(two_kap);
    mpz_ui_pow_ui(two_kap, 2, p_kap);
    
    for (int i = 0; i < p_u.size(); i++){
        mpz_init(p_y[i]);
        mpz_fdiv_q(p_y[i], p_u[i], two_kap);
        //y.push_back(p_u[i] / pow(2, p_kap));
    }
    
    mpz_clear(two_kap);
    //return y;
}

void Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    for (int i = 0; i < p_o.size(); i++){
        mpz_init_set(p_o[i], o_D.r_x[i]);
    }
    //std::vector<mpz_t> o = o_D.getDeltaList();
    //return o;
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
