//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//


// TODO
// put the y/u/z etc stuff in the correct places
// pre-generate a bunch of random ## for encoding process
// make longs??


#include "Pk.hpp"
#include "Deltas.hpp"
#include "Pri_U.hpp"
#include <iostream>

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_p(l), p_pi(1), p_q0(1), p_x0(1), p_x(tau), p_xi(l), p_ii(l), p_n(n), p_B(Theta/theta), p_s(l,std::vector<int>(Theta)), p_vert_s(Theta,std::vector<int>(l)), p_u(Theta), p_y(Theta), p_o(Theta)
{
    
    //make t state
    gmp_randstate_t p_t_state;
    gmp_randinit_mt(p_t_state);
    gmp_randseed_ui(p_t_state, time(0)); //time is not great - better TODO
    
    make_p(p_t_state);
    make_pi();
    make_q0(p_t_state);
    make_x0();
    make_x();
    make_xi();
    make_ii();
    make_s();
    make_vert_s();
    make_u();
    make_y();
    make_o();
}

//Pk::~Pk(){
    // TODO
    
    //state clear
//    gmp_randclear(p_t_state);
//}

mpz_class Pk::encode(std::vector<int> m){
    //make class state
    gmp_randclass p_class_state (gmp_randinit_mt);
    p_class_state.seed(time(0)); //TODO
    
    //m*xi
    std::vector<mpz_class> m_xi(p_l);
    for (int i = 0; i < p_l; i++){
        m_xi[i] = m[i]*p_xi[i];
    }
    
    //bi*ii
   
    std::vector<mpz_class> bi_ii(p_l);
    for (int i = 0; i < p_l; i++){
        
        mpz_class lb = power(-2,p_alphai);
        mpz_class ub = power(2,p_alphai);
        mpz_class bi = p_class_state.get_z_range(ub-lb);
        bi = bi + lb;
        bi_ii[i] = bi*p_ii[i];
    }
    
    //b*x
    std::vector<mpz_class> b_x(p_tau);
    for (int i = 0; i < p_tau; i++){

        mpz_class lb = power(-2,p_alpha);
        mpz_class ub = power(2,p_alpha);
        mpz_class b = p_class_state.get_z_range(ub-lb);
        b = b + lb;
        
        b_x[i] = b*p_x[i];
    }
    
    //summation
    mpz_class bigsum = sum_array(m_xi) + sum_array(bi_ii) + sum_array(b_x);
    
    mpz_class c = modNear(bigsum, p_x0);
    
    return c;

}

std::vector<int> Pk::decode(mpz_class c){
    std::vector<int> m(p_l);
    for (int i = 0; i < p_l; i++){
        //std::cout << "slot " << i << "\n";
        //std::cout << "pi= " << p_p[i] << "\n";
        
        mpz_class mn = modNear(c,p_p[i]);
        //std::cout << "modNear " << mn << "\n";
        
        mpz_class conv = floor_mod(mn,2);
        //std::cout << "mod2 " << conv << "\n";
        
        int i_conv = (int) conv.get_si(); //hopefully right
        m[i] = (i_conv);
    }
    return m;
}

std::vector<int> Pk::decode_squashed(mpz_class c){ //TODO gen
    std::vector<int> temp;
    return temp;
}

mpz_class Pk::recode(mpz_class c){ //TODO gen
    mpz_class temp = 0;
    return temp;
}

mpz_class Pk::H_add(mpz_class c1, mpz_class c2){
    mpz_class c = floor_mod((c1+c2),p_x0);
    return c;
}

mpz_class Pk::H_mult(mpz_class c1, mpz_class c2){
    mpz_class c = floor_mod((c1*c2),p_x0);
    return c;
}


//PRIVATE HELPER
void Pk::make_p(gmp_randstate_t p_t_state){
    for (int i = 0; i < p_l; i++){
       p_p[i] = random_prime_w(p_eta, p_t_state); //weird range 2^(n-1), 2^n
    }
}

void Pk::make_pi(){ //prod of all p[i]
    p_pi = 1;
    for (int i = 0; i < p_l; i++){
        p_pi = p_pi*p_p[i];
    }
}

void Pk::make_q0(gmp_randstate_t p_t_state){
    p_q0 = power(2,p_gam);
    mpz_class comp = floor_div(p_q0,p_pi);
    
    while (p_q0 > comp){  
        mpz_class q0_prime1 = random_prime_f0(pow(p_lam,2), p_t_state); //2^entry //need to be long?TODO
        mpz_class q0_prime2 = random_prime_f0(pow(p_lam,2), p_t_state); //2^entry
        
        p_q0 = q0_prime1*q0_prime2;
    }
}

void Pk::make_x0(){
    p_x0 = p_pi * p_q0;
}

void Pk::make_x(){ //TODO DELTAS - TODO initialize list
    Deltas x_D = Deltas(*this, p_tau, p_rhoi-1, 0);
    p_x = x_D.r_x; //getDeltaList();
}

void Pk::make_xi(){
    Deltas xi_D = Deltas(*this, p_l, p_rho, 1);
    p_xi = xi_D.r_x; //getDeltaList();
}

void Pk::make_ii(){
    Deltas ii_D = Deltas(*this, p_l, p_rho, 2);
    p_ii = ii_D.r_x; //getDeltaList();
}

void Pk::make_s(){
    for(int j = 0; j < p_l; j++){
        //all initialized to 0 originally
            p_s[j][j] = 1;
    }
    
    for(int t = 1; t < p_theta; t++){
        std::vector<int> sri = random_sample(p_B, p_l);
        for(int j = 0; j < p_l; j++){
            int k = (p_B*t)+sri[j];
            p_s[j][k] = 1;
        }
    }
    
    //print
    /*
    for (int i = 0; i < p_l; i++)
    {
        for (int j = 0; j < p_Theta; j++)
        {
            std::cout << i << ", " << j << ": " << p_s[i][j] << "\n";
        }
    }
     */
}

void Pk::make_vert_s(){
    for(int i = 0; i < p_Theta; i++){
        for(int j = 0; j < p_l; j++){
            p_vert_s[i][j] = p_s[j][i];
        }
    }
}

void Pk::make_u(){
    Pri_U priu = Pri_U(*this, p_Theta);
    p_u = priu.u_u; //.getUList();
}

void Pk::make_y(){
    mpz_class div = power(2,p_kap);
    for (int i = 0; i < p_u.size(); i++){
        p_y[i] = (p_u[i] / div); //TODO either floor div it or more likely - need to convert to rational
    }
}

void Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    p_o = o_D.r_x; //getDeltaList();
}


