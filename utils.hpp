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

mpz_class modNear(mpz_class a, mpz_class b);
mpz_class mod(mpz_class a, int b);
mpz_class mod(mpz_class a, mpz_class b);
mpz_class random_element(int l, int u);

//int set_random_seed(int seed);
std::vector<int> random_sample(int cap, int len);
//std::vector<int> sumBinary(std::vector<int> a, std::vector<int> b);
//std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
//std::vector<int> toBinary(int x, int l);
mpz_class mul_inv(mpz_class a, mpz_class b);
mpz_class CRT(std::vector<mpz_class> n, std::vector<mpz_class> a);
int kd(int i, int j);
//std::vector<int> vec_mult(int c, std::vector<int> v);
mpz_class random_prime(int l, int u);
//bool isPrime(int t);
int random_choice(std::vector<int> sample);
std::vector<int> random_sample(int range, int l);


#endif /* utils_hpp */
