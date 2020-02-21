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
#include "PseudoRandomInts.hpp"

class Pk{
private:
    //helper
    
    void make_p(gmp_randstate_t p_t_state);
    void make_pi();
    void make_q0(gmp_randstate_t p_t_state);
    void make_x0();
    

    void make_x_Delta();
    void make_xi_Delta();
    void make_ii_Delta();
    void make_o_Delta();
    void make_op_Delta();
    
    void make_s();
    void make_u();
    
    std::vector<mpz_class> make_Deltas(PseudoRandomInts r_chi, int r_lenv, int r_rho, int r_cr);
    std::vector<mpz_class> make_x_list(PseudoRandomInts chi, std::vector<mpz_class> deltas);
    std::vector<mpz_class> make_front_U(PseudoRandomInts u_pri);
    std::vector<mpz_class> make_full_u(PseudoRandomInts priu);
    std::vector<mpq_class> make_y();
    
    mpz_class recode_work(mpz_class c, std::vector<mpz_class> o);
    int permute(int j);
    std::vector<std::vector<int>> expand(mpz_class c);
    
    //Encoding Helper
    
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
    long p_seed_x;
    long p_seed_xi;
    long p_seed_ii;
    long p_seed_o;
    long p_seed_op;
    long p_seed_u;
    std::vector<mpz_class> p_x_D;
    std::vector<mpz_class> p_xi_D;
    std::vector<mpz_class> p_ii_D;
    std::vector<mpz_class> p_o_D;
    std::vector<mpz_class> p_op_D;
    int p_B;
    std::vector<std::vector<int>> p_s;
    std::vector<mpz_class> p_u_front;
    
    
    Pk(int lam, int rho, int eta, int gam, int Theta, int alpha, int tau, int l, int n=4);
    //~Pk();
    
    //TODO - DEINITIALIZE ALL THIS CRAP
    
    bool assert_parameter_correctness();
    mpz_class encode(std::vector<int> m);
    std::vector<int> decode(mpz_class c);
    std::vector<int> decode_squashed(mpz_class c);
    mpz_class recode(mpz_class c);
    mpz_class recode_and_permute(mpz_class c);
    mpz_class H_add(mpz_class c1, mpz_class c2);
    mpz_class H_mult(mpz_class c1, mpz_class c2);
    
    static Pk make_key(int size);

};


#endif /* Pk_hpp */
