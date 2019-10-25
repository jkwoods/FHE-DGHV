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

void modNear(mpz_t result, mpz_t a, mpz_t b){ //convert to long/int?
    
    //(2*a+b) / (2*b);
    mpz_t twoa;
    mpz_init_set(twoa, b);
    mpz_addmul_ui(twoa, a, 2);
    
    mpz_t twob;
    mpz_init(twob);
    mpz_mul_ui(twob, b, 2);
    
    //quotientNear
    mpz_fdiv_q(twoa, twoa, twob);
    
    //mpz_class mn = a - b*quotientNear;
    mpz_set(result, a);
    mpz_submul(result, b, twoa);
    
    mpz_clear(twoa);
    mpz_clear(twob);
}

//void mod(mpz_t result, mpz_t a, mpz_t b){
//    mpz_mod(result, a, b);
//}

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

void mul_inv(mpz_t result, mpz_t a, mpz_t b){ //TODO - finish
    mpz_t b0;
    mpz_init_set(b0, b);
    
    mpz_t x0;
    mpz_init_set_ui(x0, 0);

    mpz_t x1;
    mpz_init_set_ui(x1, 1);
    
    if (mpz_cmp_ui(b,1)==0){ //b==1
        mpz_set_ui(result, 1);
        //return 1;
    }
    mpz_t q;
    mpz_init(q);
    
    mpz_t temp_a;
    mpz_t temp_b;
    mpz_t temp_x0;
    mpz_init_set(temp_a, a);
    mpz_init_set(temp_b, b);
    
    while (mpz_cmp_ui(temp_a,1)>0){ //a>1
        mpz_fdiv_q(q, temp_a, temp_b);   //int q = a / b;
        
        mpz_set(temp_a, temp_b); //a = b;
        mpz_mod(temp_b, temp_a, temp_b); //b = a % b;
        
        mpz_set(temp_x0, x0);//int temp = x0;
        mpz_mul(x0, x0, q); //x0 = x1 - q*x0;
        mpz_sub(x0, x1, x0);
        mpz_set(x1, temp_x0); //x1 = temp;
    }
    if (mpz_cmp_ui(x1,0)<0){ //x1 < 0
        mpz_add(x1, x1, b0);
        //x1 += b0;
    }
    mpz_set(result, x1); //return x1;
    
    
    mpz_clear(b0);
    mpz_clear(x0);
    mpz_clear(x1);
    mpz_clear(q);
    mpz_clear(temp_a);
    mpz_clear(temp_b);
    mpz_clear(temp_x0);
}

void CRT(mpz_t result, std::vector<mpz_t> n, std::vector<mpz_t> a){ //chinese remainder thm
    mpz_t sum;
    mpz_init_set_ui(sum, 0);
    
    mpz_t prod;
    mpz_init_set_ui(prod, 1);
    
    for (int i = 0; i < n.size(); i++){
        mpz_mul(prod, prod, n[i]); //prod = prod*n[i];
    }
    
    mpz_t p;
    mpz_init(p);
    
    mpz_t inv;
    mpz_init(inv);
    for (int i = 0; i < n.size(); i++){
        mpz_fdiv_q(p, prod, n[i]); //mpz_class p = prod / n[i];
        
        mul_inv(inv, p, n[i]); //sum += a[i] * mul_inv(p, n[i]) * p;
        mpz_mul(inv, inv, a[i]);
        mpz_addmul(sum, inv, p);
    }
    
    mpz_mod(result, sum, prod); //return (sum % prod);
    
    mpz_clear(sum);
    mpz_clear(prod);
    mpz_clear(p);
    mpz_clear(inv);
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



