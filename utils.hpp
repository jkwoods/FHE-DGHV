//
//  utils.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <vector>
#include <gmpxx.h>
#include <math.h>
#include <gmp.h>
#include <iostream>
#include <random>


mpq_class mod_2_f(mpq_class a);
void print_vec(std::vector<int> p);
void print_vec(std::vector<mpz_class> p);
std::vector<int> to_binary(long a, int bits);
mpz_class floor_mod(mpz_class a, int b);
mpz_class floor_mod(mpz_class a, mpz_class b);
mpz_class floor_div(mpz_class a, mpz_class b);
mpz_class modNear(mpz_class a, mpz_class b);
mpz_class random_element(int l, int u);
mpz_class mul_inv(mpz_class a, mpz_class b);
mpz_class CRT(std::vector<mpz_class> n, std::vector<mpz_class> a);
int kd(int i, int j);
mpz_class power(long base, long exp);
mpz_class power(int base, int exp);
int seed_state(int seed, gmp_randclass rand_state); //TODO
mpz_class random_prime_w(int ub, gmp_randstate_t rand_state); //weird range entry
mpz_class random_prime_f0(int ub, gmp_randstate_t rand_state); //lb == 0 //2^entry
int random_choice(std::vector<int> sample);
std::vector<int> random_sample(int range, int l);
mpz_class sum_array(std::vector<mpz_class> a);
std::vector<mpz_class> sum_binary(std::vector<mpz_class> a, std::vector<mpz_class> b, mpz_class x0);

//std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);





#endif /* utils_hpp */
