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
   /*
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
    
  
    for(int i = 0; i < pk_a.p_x.size(); i++){
        if (pk_a.p_x[i] < 0){
            std::cout << "shit" << "\n";
        }
        
    }
    
    std::cout << "done" << "\n";
    
    
    std::cout << pk_a.p_pi << "\n";
    
    if(pk_a.p_eta >= (pk_a.p_alphai + pk_a.p_rhoi+1+pk_a.p_logl)){
        
        std::cout << "true" << "\n";
    }
       
    if(pk_a.p_eta >= (pk_a.p_lam * pow(log(pk_a.p_lam),2))){
           
        std::cout << "true" << "\n";
    }
    
    
    
    std::vector<int> z_vec = {1,1,1,1,1,1,1,1,1,1};
    Encoding zero = Encoding(pk_a, z_vec);
    
     std::cout << "e= " << zero.e_val << "\n";
    
    std::vector<int> z_dec = zero.decode();
    
   
    
    //print
    for (int i = 0; i < l; i++)
    {
        //std::cout << z_dec[i] << "\n";
    }
    
    
    // timing tests
    std::cout << "Tests\n";
    std::cout << time(0);
    
    std::cout << "End\n";
    
*/

    mpz_class a = -8;
    mpz_class b = -1;
    mpz_class c = -47;
    mpz_class d = 0;
    
    mpz_class t1 = a % 2;
    mpz_class t2 = b % 2;
    mpz_class t3 = c % 2;
    mpz_class t4 = d % 2;

    
    std::cout << t1 << "\n";
    std::cout << t2 << "\n";
    std::cout << t3 << "\n";
    std::cout << t4 << "\n";

    //mpz_class t1 = modNear(a, b);
    //mpz_class t2 = modNear(b, a);
    
    //mpz_class t1 = mul_inv(a,b);
    //mpz_class t2 = mul_inv(b,a);
    
    
    return 0;
}

