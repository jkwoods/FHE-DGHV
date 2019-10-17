//
//  utils.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <stdio.h>
#include <vector>

int modNear(int a, int b);
int mod(int a,int b);
int random_element(int l, int u);
int set_random_seed(int seed);
std::vector<int> random_sample(int cap, int len);
std::vector<int> sumBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> toBinary(int x, int l);
int mul_inv(int a, int b);
int CRT(std::vector<int> n, std::vector<int> a);
int kd(int i, int j);
std::vector<int> vec_mult(int c, std::vector<int> v);
int random_prime(int l, int u);
bool isPrime(int t);
int random_choice(std::vector<int> sample);


#endif /* utils_hpp */
