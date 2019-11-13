//
//  utils.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "utils.hpp"
#include <math.h>

mpz_class modNear(mpz_class a, mpz_class b){ //convert to long/int?
    mpz_class quotientNear = (2*a+b) / (2*b); // autmatically rounds towards floor (check?)
    mpz_class mn = a - b*quotientNear;
    return mn;
}

//void mod(mpz_t result, mpz_t a, mpz_t b){
//    mpz_mod(result, a, b);
//}

//int random_element(int l, int u){ //not power of two, lower bound assumed zero
//   int r = 1;
//    if (l >= 0){
//        r = (rand() % u) + l;
//    } else {
//        r = (rand() % ((-1*l)+u)) - l;
//   }
//    return r;
//}

/*void random_element_pow2(mpz_t elt, int l, int u, gmp_randstate_t state){
    mpz_urandomb(elt, state, l+u);
    
    mpz_t two_l;
    mpz_init(two_l);
    mpz_ui_pow_ui(two_l, 2, l);
    
    mpz_sub(elt, elt, two_l);
} */

/*
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
} */

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

//std::vector<int> xorBinary(std::vector<int> a, std::vector<int> b);
//std::vector<int> toBinary(int x, int l);

mpz_class mul_inv(mpz_class a, mpz_class b){ //TODO - finish
    mpz_class b0 = b;
    mpz_class x0 = 0;
    mpz_class x1 = 1;
    if (b == 1){
        return 1;
    }
    while (a>1){
        mpz_class q = a / b; //floor
        mpz_class temp = b;
        b = a % b;
        a = temp;
        
        temp = x0;
        x0 = x1 - q * x0;
        x1 = temp;
    }
    if (x1 < 0){
        x1 += b0;
    }
    return x1;
}

mpz_class CRT(std::vector<mpz_class> n, std::vector<mpz_class> a){ //chinese remainder thm
    mpz_class prod = 1;
    for (int i = 0; i < n.size(); i++){
        prod = prod*n[i];
    }
    
    mpz_class sum = 0;
    for (int i = 0; i < n.size(); i++){
        mpz_class p = prod / n[i]; //floor
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

//std::vector<int> vec_mult(int c, std::vector<int> v){
//  std::vector<int> r;
//  for (int i = 0; i < v.size(); i++) {
//    r.push_back(c*v.at(i));
//  }
//  return r;
//}

//int random_prime(int l, int u){
//    int r = random_element(l, u);
//    while (!isPrime(r)){
//        r = random_element(l, u);
//    }
//    return r;
//}

//bool isPrime(int t){
//    for(int i = 2; i < sqrt(t); i++){
//       if (t % i == 0){
//            return false;
//        }
//    }
//    return true;
//}

int random_choice(std::vector<int> sample){
    int r = random_element(0, sample.size()-1); //TODO need -1??
    return sample[r];
}

std::vector<int> random_sample(int range, int l){
    std::vector<int> sample;
    for(int i = 0; i < range; i++){
        sample.push_back(i);
    }
    std::random_shuffle(sample.begin(), sample.end());
    std::vector<int> cut_sample;
    for(int i = 0; i < l; i++){
        cut_sample.push_back(sample[i]);
    }
    return cut_sample;
}

