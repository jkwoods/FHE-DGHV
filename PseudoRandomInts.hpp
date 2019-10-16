//
//  PseudoRandomInts.hpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef PseudoRandomInts_hpp
#define PseudoRandomInts_hpp

#include <stdio.h>
#include <vector>

class PseudoRandomInts{
public:
    PseudoRandomInts(int x0, int len);
    int r_len;
    
    int getSeed();
    std::vector<int> getList();
    std::vector<int> makeList();

private:
    int r_x0;
    int r_seed;
    std::vector<int> r_list;
    
};

#endif /* PseudoRandomInts_hpp */
