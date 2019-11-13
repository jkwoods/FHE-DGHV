//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//


// TODO - declare vector size at beginning, change rand (make state, etc), change accumulate, change pow
// put the y/u/z etc stuff in the correct places


#include "Pk.hpp"
#include "Deltas.hpp"
#include "Pri_U.hpp"

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_p(l), p_pi(), p_q0(), p_x0(), p_x(tau), p_xi(l), p_ii(l), p_n(n), p_B(Theta/theta), p_s(p_l, std::vector<int> (p_Theta)), p_vert_s(p_Theta, std::vector<int> (p_l)), p_u(), p_y(), p_o(Theta) {
    
    make_p();
    make_pi();
    make_q0();
    make_x0();
    make_x();
    make_xi();
    make_ii();
    make_s();
    make_vert_s();
    make_u();
    make_y();
    make_o();
    
    make_state();
}

Pk::~Pk(){
    // TODO
}

mpz_class Pk::encode(std::vector<int> m){
    //m*xi
    std::vector<mpz_class> m_xi(p_l);
    for (int i = 0; i < p_l; i++){
        m_xi.push_back(m[i]*p_xi[i]);
    }
    
    //bi*ii
   
    std::vector<mpz_class> bi_ii(p_l);
    for (int i = 0; i < p_l; i++){
        
        random_element_pow2(bi, p_alphai, p_alphai, p_state); //needs to raise these to -2 and 2 //TODO - call state before all rand numbers
        
        mpz_class bi = (random_element(pow(-2,p_alphai),pow(2,p_alphai)));
        bi_ii.push_back(bi*p_ii[i]);
    }
    
    //b*x
    std::vector<mpz_class> b_x(p_tau);
    for (int i = 0; i < p_tau; i++){
        random_element_pow2(b, p_alpha, p_alpha, p_state); //needs to raise these to -2 and 2
        
        mpz_class b = (random_element(pow(-2,p_alpha),pow(2,p_alpha)));
        b_x.push_back(b*p_x[i]);
    }
    
    //summation
    mpz_class bigsum = accumulate(m_xi.begin(), m_xi.end(), 0) + accumulate(bi_ii.begin(), bi_ii.end(), 0) + accumulate(b_x.begin(), b_x.end(), 0);
    
    mpz_class c = modNear(bigsum, p_x0);
    
    return c;

}

std::vector<int> Pk::decode(mpz_class c){
    std::vector<int> m(p_l);
    for (int i = 0; i < p_l; i++){
        mpz_class conv = mod(modNear(c,p_p[i]),2);
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
}

mpz_class Pk::H_add(mpz_class c1, mpz_class c2){
    mpz_class c = mod(c1+c2,p_x0);
    return c;
}

mpz_class Pk::H_mult(mpz_class c1, mpz_class c2){
    mpz_class c = mod(c1*c2,p_x0);
    return c;
}

//private helper
void Pk::make_p(){
    for (int i = 0; i < p_l; i++){
       p_p[i] = (random_prime(pow(2,p_eta-1), pow(2,p_eta)));
    }
}

void Pk::make_pi(){ //prod of all p[i]
    p_pi = 1;
    for (int i = 0; i < p_l; i++){
        p_pi = p_pi*p_p[i];
    }
}

void Pk::make_q0(){
    p_q0 = pow(2,p_gam); //TODO
    mpz_class comp = p_q0 / p_pi;
    
    while (p_q0 > comp){     //while (q0 > (pow(2,p_gam)/p_pi)){
        int q0_prime1 = random_prime(0, pow(2,pow(p_lam,2)));
        int q0_prime2 = random_prime(0, pow(2,pow(p_lam,2)));
        
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
    for(int i = 0; i < p_l; i++){
        for (int j = 0; j < p_Theta; j++){
            std::vector<int> fill(p_B, 0); //initialize all to 0
            if (j==0){
                fill[j] = 1;
                p_s[i].insert(std::end(s[i]), std::begin(fill), std::end(fill));
            } else {
                p_s[i].insert(std::end(s[i]), std::begin(fill), std::end(fill));
            }
        }
    }
    
    for(int i = 0; i < p_theta; i++){ //TODO??
        std::vector<int> sri = random_sample(p_B, p_l);
        for(int j = 0; i < p_l; i++){
            int k = (p_B*i)+sri[j];
            p_s[j][k] = 1;
        }
    }
}

void Pk::make_vert_s(){
    for(int i = 0; i < p_Theta; i++){
        for(int j = 0; j < p_l; j++){
            p_vert_s[i][j] = p_s[j][i];
        }
    }
}

void Pk::make_u(){
    Pri_U priu = Pri_U(*this);
    p_u = priu.u_u; //.getUList();
}

void Pk::make_y(){
    mpz_class div = pow(2, p_kap); //TODO
    for (int i = 0; i < p_u.size(); i++){
        p_y[i] = (p_u[i] / div);
    }
}

void Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    p_o = o_D.r_x; //getDeltaList();
}

void Pk::make_state(){ //TODO
    gmp_randinit_mt(p_state);
    gmp_randseed_ui(p_state, time(0)); //time as seed TODO - check this shit
}

