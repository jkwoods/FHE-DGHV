//
//  main.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//

#include <iostream>
#include <vector>
#include <bitset>
#include "Pk.hpp"
#include "Encoding.hpp"
#include "utils.hpp"
#include <ctime>

Encoding encode(std::vector<int> m, Pk pk);


int main(int argc, const char * argv[]) {
    
    // make pk
    std::cout << "Making Keys\n";
    int lam=12;
    int rho=26;
    int rhoi=26;
    int eta=1988;
    int gam=147456;
    int Theta=150;
    int theta=15;
    int kap=149446;
    int alpha=936;
    int alphai=936;
    int tau=188;
    int l=10;
    
    Pk pk_a = Pk(lam, rho, rhoi, eta, gam, Theta, theta, kap, alpha, alphai, tau, l);
    
    //TODO first - check all utils, go from there
    
    std::vector<int> z_vec = {0,0,0,0,0,0,0,0,0,0};
    Encoding zero = Encoding(pk_a, z_vec);
    
    std::vector<int> z_dec = zero.decode();
    
    //print
    for (int i = 0; i < l; i++)
    {
        std::cout << z_dec[i] << "\n";
    }
    
    
    // timing tests
    std::cout << "Tests\n";
    std::cout << time(0);
    
    std::cout << "End\n";
    
    return 0;
}

