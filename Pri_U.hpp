//
//  Pri_U.hpp
//  FHE
//
//  Created by Woods, Jess on 10/17/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef Pri_U_hpp
#define Pri_U_hpp

#include <stdio.h>
#include <vector>
#include "Pk.hpp"
#include "PseudoRandomInts.hpp"

class Pri_U{
public:
    Pri_U(Pk pk);
    
    std::vector<int> getUList();
    PseudoRandomInts u_pri;
    
private:
    std::vector<int> u_u;
    Pk u_pk;
    
    std::vector<int> makeU();
    PseudoRandomInts makePri();
};


#endif /* Pri_U_hpp */
