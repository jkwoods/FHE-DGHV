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
#include <numeric>
#include <gmp.h>
#include "Encoding.hpp"

class Pk{
private:
    //helper
    void make_p();
    void make_pi();
    void make_q0();
    void make_x0();
    void make_x();
    void make_xi();
    void make_ii();
    std::vector<std::vector<int>> make_s();
    std::vector<std::vector<int>> make_vert_s();
    void make_u();
    void make_y();
    void make_o();
    std::vector<int> random_sample(int range, int l);
    
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
    std::vector<mpz_t> p_p;
    mpz_t p_pi;
    mpz_t p_q0;
    mpz_t p_x0;
    std::vector<mpz_t> p_x;
    std::vector<mpz_t> p_xi;
    std::vector<mpz_t> p_ii;
    int p_B;
    std::vector<std::vector<int>> p_s;
    std::vector<std::vector<int>> p_vert_s;
    std::vector<mpz_t> p_u;
    std::vector<mpz_t> p_y;
    std::vector<mpz_t> p_o;
    
    
    Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n=4);
    
    //TODO - DEINITIALIZE ALL THIS CRAP
    
    void encode(mpz_t c, std::vector<int> m);
    std::vector<int> decode(mpz_t c);
    std::vector<int> decode_squashed(mpz_t c);
    void recode(mpz_t r, mpz_t c);
    void H_add(mpz_t added, mpz_t c1, mpz_t c2);
    void H_mult(mpz_t multed, mpz_t c1, mpz_t c2);
    

};


#endif /* Pk_hpp */
