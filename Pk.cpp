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
#include "utils.hpp"

Pk::Pk(int lam, int rho, int eta, int gam, int Theta, int alpha, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rho+lam), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(Theta/l), p_kap(64*(gam/64+1)-1), p_alpha(alpha), p_alphai(alpha+lam), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_p(l), p_pi(1), p_q0(1), p_x0(1), p_x(tau), p_xi(l), p_ii(l), p_n(n), p_B(l), p_s(l,std::vector<int>(Theta)), p_vert_s(Theta,std::vector<int>(l)), p_u(Theta), p_y(Theta), p_o(Theta) //for kap, c++ trucates (rounds down) so its all gud
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

std::vector<std::vector<int>> Pk::expand(mpz_class c){
    //TODO - recover u
    std::vector<std::vector<int>> z(p_Theta,std::vector<int>(p_n+1));
    for (int i = 0; i < p_Theta; i++){
        mpq_class zi = c * p_y[i];
        mpq_class mod = mod_2_f(zi);
        //std::cout << mpf_class(mod) << "\n";
        
        //convert each zi to vector of binary
        //mult by 2^n (bits of precision)
        mpq_class zi_mult = mod * (pow(2,p_n));
        //std::cout << mpf_class(zi_mult) << "\n";
        
        //convert to int (cut off ends)
        mpz_class zi_conv = mpz_class(zi_mult);
        if ((zi_mult - zi_conv) > 0.5){
            zi_conv += 1;
        }
        //std::cout << zi_conv << "\n";
        
        if (zi_conv.fits_sint_p() == 0){ //doesn't fit
            std::cout << "Error in generation of z\n";
        }
        int zi_round = zi_conv.get_si(); //round
        //std::cout << zi_round << "\n";
        
        //put in binary
        std::vector<int> ex_z(p_n+1);
        ex_z = to_binary(zi_round, p_n+1); //[lsb, ...., msb] max 31 (can we ever possibly get a 32? what then) TODO
        //print_vec(ex_z);
        
        z[i] = ex_z;
    }
    return z;
}

mpz_class Pk::recode(mpz_class c){
    std::vector<std::vector<int>> z = expand(c);
    
    std::vector<std::vector<mpz_class>> z_mult(p_Theta,std::vector<mpz_class>(p_n+1));
    
    //z * sk - correct
    for(int i = 0; i < p_Theta; i++){
        for(int j = 0; j < p_n+1; j++){
            z_mult[i][j] = (z[i][j] * p_o[i]);
        }
    }
    
    //add to make as
    std::vector<mpz_class> a(p_n+1);
    for(int i = 0; i < p_Theta; i++){
        a = sum_binary(a,z_mult[i]);

        for (int j = 0; j < p_n+1; j++){
            //a[j] = floor_mod(a[j], p_x0);
            //std::cout << "a; i= " << i << "\n";
            //print_vec(decode(a[j]));
        }
        
    }
    
    //print_vec(decode(a[a.size()-2]));
    //std::cout << "c & 1= " << (c & 1) << "\n";
    
    mpz_class two = a[a.size()-1] + a[a.size()-2]; //correct??
    mpz_class c_prime = two + (c & 1);
    
    return c_prime;
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
    for (int j = 0; j < 15; j++)
    {
        for (int i = 0; i < p_l; i++)
        {
            std::cout << p_s[i][(j*10)+0] << p_s[i][(j*10)+1] << p_s[i][(j*10)+2] << p_s[i][(j*10)+3] << p_s[i][(j*10)+4] << p_s[i][(j*10)+5] << p_s[i][(j*10)+6] << p_s[i][(j*10)+7] << p_s[i][(j*10)+8] << p_s[i][(j*10)+9] << "\n";
        }
        std::cout << "\n\n";
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
        p_y[i] = mpq_class(p_u[i],div); //rational u[i]/(2^kap)
        p_y[i].canonicalize();
    }
}

void Pk::make_o(){
    Deltas o_D = Deltas(*this, p_Theta, p_rho, 3);
    p_o = o_D.r_x; //getDeltaList();
}

bool Pk::assert_parameter_correctness(){
    bool a = p_rho >= 2*p_lam; //brute force noise attack
    std::cout << a << "\n";
    bool b = p_eta >= p_alphai + p_rhoi + 1 + log2(p_l); // correct decoding
    std::cout << b << "\n";
    bool c = p_eta >= p_rho * (p_lam*(pow(log(p_lam),2))); //squashed decode circut
    std::cout << c << "\n";
    bool d = p_gam > pow(p_eta, 2) * log(p_lam); //lattice attack
    std::cout << d << "\n";
    bool e = (p_alpha * p_tau) >= p_gam + p_lam; //leftover hash lemma
    std::cout << e << "\n";
    bool f = p_tau >= p_l * (p_rhoi + 2) + p_lam; //leftover hash lemma
    std::cout << f << "\n";
    bool g = (p_Theta % p_l == 0);
    
    return a && b && c && d && e && f && g;
}

Pk Pk::make_key(int size){
    if (size == 0){
        int lam=42;
        int rho=84;
        int eta=1940;
        int gam=140400;
        int Theta=150;
        int alpha=1764;
        int tau=74088;
        int l=10;
        
        return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
    } else if (size == 1){
        int lam=52;
        int rho=104;
        int eta=146767;
        int gam=21550;
        int Theta=555;
        int alpha=1476;
        int tau=5900;
        int l=37;
        
        return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
    } else if (size == 2){
        int lam=52;
        int rho=104;
        int eta=146767;
        int gam=21550;
        int Theta=555;
        int alpha=1476;
        int tau=5900;
        int l=37;
        
        return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
    } else if (size == 3){
        int lam=52;
        int rho=104;
        int eta=146767;
        int gam=21550;
        int Theta=555;
        int alpha=1476;
        int tau=5900;
        int l=37;
        
        return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
    } else {
        std::cout << "Size not correctly specified. Small key being made.\n";
        int lam=52;
        int rho=104;
        int eta=146767;
        int gam=21550;
        int Theta=555;
        int alpha=1476;
        int tau=5900;
        int l=37;
        
        return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
        
    }
}

