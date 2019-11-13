//
//  Pk.hpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//

#ifndef Pk_hpp
#define Pk_hpp

#include <stdio.h>
#include <vector>
#include "utils.hpp"
#include <math.h>
#include <cmath>
#include <gmpxx.h>

class Pk{
private:
    //helper
    
    void make_p(gmp_randstate_t p_t_state);
    void make_pi();
    void make_q0(gmp_randstate_t p_t_state);
    void make_x0();
    void make_x();
    void make_xi();
    void make_ii();
    void make_s();
    void make_vert_s();
    void make_u();
    void make_y();
    void make_o();
    
public:
    int p_lam;
    int p_rho;
    int p_rhoi;
    int p_eta;
    int p_gam;
    int p_Theta;
    int p_theta;
    int p_n;
    int p_kap;
    int p_alpha;
    int p_alphai;
    int p_tau;
    int p_l;
    int p_logl;
    std::vector<mpz_class> p_p;
    mpz_class p_pi;
    mpz_class p_q0;
    mpz_class p_x0;
    std::vector<mpz_class> p_x;
    std::vector<mpz_class> p_xi;
    std::vector<mpz_class> p_ii;
    int p_B;
    std::vector<std::vector<int>> p_s;
    std::vector<std::vector<int>> p_vert_s;
    std::vector<mpz_class> p_u;
    std::vector<mpz_class> p_y;
    std::vector<mpz_class> p_o;
    
    
    Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n=4);
    //~Pk();
    
    //TODO - DEINITIALIZE ALL THIS CRAP
    
    mpz_class encode(std::vector<int> m);
    std::vector<int> decode(mpz_class c);
    std::vector<int> decode_squashed(mpz_class c);
    mpz_class recode(mpz_class c);
    mpz_class H_add(mpz_class c1, mpz_class c2);
    mpz_class H_mult(mpz_class c1, mpz_class c2);
    

};


#endif /* Pk_hpp */
