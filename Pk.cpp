//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//


// TODO - declare vector size at beginning, change rand (make state, etc), change accumulate
// put the y/u/z etc stuff in the correct places


#include "Pk.hpp"
#include "Deltas.hpp"
#include "Pri_U.hpp"

Pk::Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rhoi), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(theta), p_kap(kap), p_alpha(alpha), p_alphai(alphai), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_p(make_p()), p_pi(make_pi()), p_q0(make_q0()), p_x0(make_x0()), p_x(make_x()), p_xi(make_xi()), p_ii(make_ii()), p_n(n), p_B(Theta/theta), p_s(make_s()), p_vert_s(make_vert_s()), p_u(make_u()), p_y(make_y()), p_o(make_o()) {
    

    make_state();
}

Pk::~Pk(){
    // TODO
}

mpz_class Pk::encode(std::vector<int> m){
    //m*xi
    std::vector<mpz_class> m_xi;
    for (int i = 0; i < p_l; i++){
        m_xi.push_back(m[i]*p_xi[i]);
    }
    
    //bi*ii
   
    std::vector<mpz_class> bi_ii;
    for (int i = 0; i < p_l; i++){
        
        random_element_pow2(bi, p_alphai, p_alphai, p_state); //needs to raise these to -2 and 2 //TODO - call state before all rand numbers
        
        mpz_class bi = (random_element(pow(-2,p_alphai),pow(2,p_alphai)));
        bi_ii.push_back(bi*p_ii[i]);
    }
    
    //b*x
    std::vector<mpz_class> b_x;
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
std::vector<mpz_class> Pk::make_p(){
    std::vector<mpz_class> p;
    for (int i = 0; i < p_l; i++){
       p.push_back(random_prime(pow(2,p_eta-1), pow(2,p_eta)));
    }
    return p;
}

mpz_class Pk::make_pi(){ //prod of all p[i]
    mpz_class pi = 1;
    for (int i = 0; i < p_l; i++){
        pi = pi*p_p[i];
    }
    return pi;
}

mpz_class Pk::make_q0(){
    mpz_class q0 = pow(2,p_gam); //TODO
    mpz_class comp = q0 / p_pi;
    
    while (q0 > comp){     //while (q0 > (pow(2,p_gam)/p_pi)){
        int q0_prime1 = random_prime(0, pow(2,pow(p_lam,2)));
        int q0_prime2 = random_prime(0, pow(2,pow(p_lam,2)));
        
        q0 = q0_prime1*q0_prime2;
    }
    
    return q0;
}

mpz_class Pk::make_x0(){
    mpz_class x0 = p_pi * p_q0;
    return x0;
}

std::vector<mpz_class> Pk::make_x(){ //TODO DELTAS - TODO initialize list
    Deltas x_D = Deltas(*this, p_tau, p_rhoi-1, 0);
    std::vector<mpz_class> x = x_D.r_x; //getDeltaList();
    return x;
}

std::vector<mpz_class> Pk::make_xi(){
    Deltas xi_D = Deltas(*this, p_l, p_rho, 1);
    std::vector<mpz_class> xi = xi_D.r_x; //getDeltaList();
    return xi;
}

std::vector<mpz_class> Pk::make_ii(){
    Deltas ii_D = Deltas(*this, p_l, p_rho, 2);
    std::vector<mpz_class> ii = ii_D.r_x; //getDeltaList();
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

std::vector<mpz_class> Pk::make_u(){
    Pri_U priu = Pri_U(*this);
    std::vector<mpz_class> u = priu.u_u; //.getUList();
    return u;
}

std::vector<mpz_class> Pk::make_y(){
    std::vector<mpz_class> y;
    mpz_class div = pow(2, p_kap); //TODO
    for (int i = 0; i < p_u.size(); i++){
        y.push_back(p_u[i] / div);
    }
    return y;
}

std::vector<mpz_class> Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    std::vector<mpz_class> o = o_D.r_x; //getDeltaList();
    return o;
}

void Pk::make_state(){ //TODO
    gmp_randinit_mt(p_state);
    gmp_randseed_ui(p_state, time(0)); //time as seed TODO - check this shit
}

