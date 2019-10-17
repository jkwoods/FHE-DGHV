//
//  RandomPrimeGenerator.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "RandomPrimeGenerator.hpp"
#include "math.h"

RandomPrimeGenerator::RandomPrimeGenerator(int lower, int upper): p_upper(upper), p_lower(lower){}

int RandomPrimeGenerator::getPrime(){
    int r = rand();
    while (!isPrime(r)){
        r = rand();
    }
}

bool isPrime(int t){
    for(int i = 2; i < sqrt(t); i++){
        if (t % i == 0){
            return false;
        }
    }
    return true;
}

