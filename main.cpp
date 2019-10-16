//
//  main.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include <iostream>
#include <vector>
#include <bitset>
#include "Pk.hpp"
#include "Encoding.hpp"
#include "utils.hpp"
#include "RandomPrimeGenerator.hpp"
#include <ctime>

Encoding encode(std::vector<int> m, Pk pk);


int main(int argc, const char * argv[]) {
    
    // make pk
    std::cout << "Making Keys\n";
    
    
    // timing tests
    std::cout << "Tests\n";
    std::cout << time(0);
    
    
    
    
    std::cout << "End\n";
    
    return 0;
}

Encoding encode(std::vector<int> m, Pk pk){
    int c = pk.encode(m);
    Encoding ec = Encoding(c,pk);
    return ec;
}

