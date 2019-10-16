//
//  RandomPrimeGenerator.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef RandomPrimeGenerator_hpp
#define RandomPrimeGenerator_hpp

#include <stdio.h>
#include <vector>

class RandomPrimeGenerator{
public:
    RandomPrimeGenerator(int limit);
    int getPrime();

private:
    std::vector<int> p_list;
    std::vector<int> p_limit;
    int p_index;
    
    std::vector<int> makeList();

};

#endif /* RandomPrimeGenerator_hpp */
