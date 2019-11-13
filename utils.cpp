//
//  utils.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "utils.hpp"
#include <math.h>
#include <gmp.h>
#include <iostream>

mpz_class modNear(mpz_class a, mpz_class b){ //convert to long/int?
    mpz_class quotientNear = (2*a+b) / (2*b); // autmatically rounds towards floor (check?)
    mpz_class mn = a - b*quotientNear;
    return mn;
}

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

mpz_class power(int base, int exp){
    mpz_t p;
    mpz_init(p);
    if (exp < 0){
         std::cout << "something is going wrong - negative exp in power function";
    }
    if (base < 0){ //handling for negative base (-2 mostly)
        base = base * (-1);
        mpz_ui_pow_ui(p, base, exp);

        mpz_class p_class(p);
        mpz_clear(p);
        
        return p_class*(-1);
    } else {
        mpz_ui_pow_ui(p, base, exp);

        mpz_class p_class(p);
        mpz_clear(p);
        
        return p_class;
    }

}

mpz_class random_prime_w(int ub, gmp_randstate_t rand_state){ //wierd range
    //generate mpz_t rand
    mpz_t p;
    mpz_init(p);
    mpz_urandomb(p, rand_state, ub-1); //0 - 2^(n-1)
    
    mpz_t ub_pow;
    mpz_init(ub_pow);
    mpz_ui_pow_ui(ub_pow, 2, ub-1);
    
    mpz_add(p, p, ub_pow);// + 2^(n-1)
    
    //check if prime
    while(mpz_probab_prime_p(p, 30) == 0){ //not prime
        mpz_urandomb(p, rand_state, ub-1); //0 - 2^(n-1)
        mpz_add(p, p, ub_pow);// + 2^(n-1)
    }
    
    //cast as mpz_class
    mpz_class p_class(p);
    
    mpz_clear(p);
    mpz_clear(ub_pow);
    
    return p_class;
}

mpz_class random_prime_f0(int ub, gmp_randstate_t rand_state){ //2^entry (only used one place)
    //generate mpz_t rand
    mpz_t p;
    mpz_init(p);
    mpz_urandomb(p, rand_state, ub);
    
    //check if prime
    while(mpz_probab_prime_p(p, 30) == 0){ //not prime
        mpz_urandomb(p, rand_state, ub);
    }
    
    //cast as mpz_class
    mpz_class p_class(p);
    mpz_clear(p);
    return p_class;
}


/*
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
 */


int random_choice(std::vector<int> sample){ //TODO test
    srand(time(0));
    int r = rand() % sample.size();
    return sample[r];
}

std::vector<int> random_sample(int range, int l){
    std::vector<int> sample(range);
    for(int i = 0; i < range; i++){
        sample[i] = i;
    }
    std::random_shuffle(sample.begin(), sample.end());
    std::vector<int> cut_sample(l);
    for(int i = 0; i < l; i++){
        cut_sample[i] = sample[i];
    }
    return cut_sample;
}

mpz_class sum_array(std::vector<mpz_class> a){
    mpz_class suma = 0;
    for(int i = 0; i < a.size(); i++){
        suma = suma+a[i];
    }
    return suma;
    
    
}
