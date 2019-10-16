//
//  RandomPrimeGenerator.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "RandomPrimeGenerator.hpp"

RandomPrimeGenerator::RandomPrimeGenerator(int limit): p_list(makeList()), p_limit(limit), p_index(0){}

int RandomPrimeGenerator::getPrime(){
    int r = p_list[p_index];
    p_index++;
    
    return r;
}

std::vector<int> RandomPrimeGenerator::makeList(){
    std::vector<int> fake;
    for(int i = 0; i < 20; i++){
        srand(time(0));
        int r = rand();
    //    while(){
            //TODO FUCK
    //    }
        
        
    }
    
    return fake;
}
