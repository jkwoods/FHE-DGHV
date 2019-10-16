//
//  utils.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "utils.hpp"

int modNear(int a, int b){
  int quotientNear = (2*a+b) / (2*b);
  int mn = a - b*quotientNear;
  return mn;
}

int mod(int a, int b){
  if (b == 2){
    return (a & 1);
  } else {
    return (a % b);
  }
}

int random_prime(int l, int u);
int random_element(int l, int u);
int set_random_seed(int seed);
std::vector<int> sumBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> toBinary(int x, int l);
int mul_inv(int a, int b);
int CRT(int n, int a);

int kd(int i, int j){
  if (i == j){
    return 1;
  } else {
    return 0;
  }
}

std::vector<int> vec_mult(int c, std::vector<int> v){
  std::vector<int> r;
  for (int i = 0; i < v.size(); i++) {
    r.push_back(c*v.at(i));
  }
  return r;
}


