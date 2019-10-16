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

Encoding encode(std::vector<int> m, Pk pk);


int main(int argc, const char * argv[]) {
    // insert code here...
    
    
    
    
    
    std::cout << "Tests\n";
    
    return 0;
}

Encoding encode(std::vector<int> m, Pk pk){
    int c = pk.encode(m);
    Encoding ec = Encoding(c,pk);
    return ec;
}

