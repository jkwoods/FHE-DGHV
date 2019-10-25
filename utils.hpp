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
#include <gmp.h>

void modNear(mpz_t result, mpz_t a, mpz_t b);
//int mod(int a,int b);
int random_element(int l, int u);
void random_element_pow2(mpz_t elt, int l, int u);
int set_random_seed(int seed);
std::vector<int> random_sample(int cap, int len);
std::vector<int> sumBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> toBinary(int x, int l);
int mul_inv(mpz_t result, mpz_t a, mpz_t b);
int CRT(mpz_t result, std::vector<mpz_t> n, std::vector<mpz_t> a);
int kd(int i, int j);
std::vector<int> vec_mult(int c, std::vector<int> v);
int random_prime(int l, int u);
bool isPrime(int t);
int random_choice(std::vector<int> sample);
std::vector<int> random_sample(int range, int l);


#endif /* utils_hpp */
