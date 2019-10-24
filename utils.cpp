//
//  utils.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "utils.hpp"
#include <math.h>
#include <stdio.h>
#include <gmp.h>

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

int random_element(int l, int u){
    int r = 1;
    if (l >= 0){
        r = (rand() % u) + l;
    } else {
        r = (rand() % ((-1*l)+u)) - l;
    }
    return r;
}

int set_random_seed(int seed){ //if seed = 0, randomize and return it, else use seed
    int s = seed;
    if (seed == 0){
        srand(time(0));
        s = rand();
        srand(s);
    } else {
        srand(s);
    }
    return s;
}

std::vector<int> sumBinary(std::vector<int> a, std::vector<int> b){
    std::vector<int> c;
    c.push_back(a[0]+b[0]);
    int carry = a[0]*b[0];
    for (int i = 1; i < a.size()-1; i++){
        int carry2 = (a[i]+b[i])*carry+a[i]*b[i];
        c.push_back(a[i]+b[i]+carry);
        carry = carry2;
    }
    c.push_back(a[-1]+b[-1]+carry);
    return c;
}

std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
std::vector<int> toBinary(int x, int l);

int mul_inv(int a, int b){
    int b0 = b;
    int x0 = 0;
    int x1 = 1;
    if (b==1){
        return 1;
    }
    while (a>1){
        int q = a / b;
        a = b;
        b = a % b;
        int temp = x0;
        x0 = x1 - q*x0;
        x1 = temp;
    }
    if (x1 < 0){
        x1 += b0;
    }
    return x1;
}

int CRT(std::vector<int> n, std::vector<int> a){ //chinese remainder thm
    int sum = 0;
    int prod = 1;
    for (int i = 0; i < n.size(); i++){
        prod = prod*n[i];
    }
    for (int i = 0; i < n.size(); i++){
        int p = prod / n[i];
        sum += a[i] * mul_inv(p, n[i]) * p;
    }
    return (sum % prod);
}

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

int random_prime(int l, int u){
    int r = random_element(l, u);
    while (!isPrime(r)){
        r = random_element(l, u);
    }
    return r;
}

bool isPrime(int t){
    for(int i = 2; i < sqrt(t); i++){
        if (t % i == 0){
            return false;
        }
    }
    return true;
}

int random_choice(std::vector<int> sample){
    int r = random_element(0, sample.size()-1); //TODO need -1??
    return sample[r];
}



