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
#include "MultEncodings.hpp"
#include <ctime>
#include <omp.h>

/* TODO
- Optimize recode
 

 */

int main(int argc, const char * argv[]) {

    
     int lam=12;
     int rho=26;
     int eta=1988;
     int gam=147456;
     int Theta=150;
     int alpha=936;
     int tau=188;
     int l=10;


     std::cout << "Making Keys\n";
     Pk pk_a = Pk(lam, rho, eta, gam, Theta, alpha, tau, l);

    //bool test = pk_a.assert_parameter_correctness();
    //std::cout << test << "\n";
    
    
    std::vector<int> z_vec = {0,0,0,0,0,0,0,0,0,0}; //,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


    double start = omp_get_wtime();
    Encoding zero = Encoding(pk_a, z_vec);
    double end = omp_get_wtime();    
    std::cout << "encoding time: " << (end-start) << "\n";

    std::vector<int> o_vec = {1,1,1,1,1,1,1,1,1,1}; //1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    Encoding one = Encoding(pk_a, o_vec);
    
    //std::vector<int> a_vec = {1,1,1,1,1,0,0,0,0,0};
    //Encoding aa = Encoding::selector(a_vec,zero,one); //Encoding(pk_a, a_vec);
    
    //std::vector<int> b_vec = {0,1,0,1,0,1,0,1,0,1};
    //Encoding bb = zero.neg(); //Encoding(pk_a, b_vec);
    
    
    Encoding a = one+zero;
    Encoding b = one*zero;
    
    start = omp_get_wtime();
    std::vector<int> z_dec = a.decode();
    end = omp_get_wtime():
    std::cout << "decoding time: " << (end-start) << "\n";
    std::vector<int> o_dec = b.decode();
    
    Encoding bad = b*zero;
    std::vector<int> a_dec = bad.decode();
    
    start = omp_get_wtime();
    b.recode();
    end = omp_get_wtime():
    std::cout << "recoding time: " << (end-start) << "\n";

    Encoding good = b*zero;
    std::vector<int> b_dec = good.decode();
    
    //print
    
    std::cout << z_dec[0] << z_dec[1] << z_dec[2] << z_dec[3] << z_dec[4] << z_dec[5] << z_dec[6] << z_dec[7] << z_dec[8] << z_dec[9] <<  "\n";
    std::cout << o_dec[0] << o_dec[1] << o_dec[2] << o_dec[3] << o_dec[4] << o_dec[5] << o_dec[6] << o_dec[7] << o_dec[8] << o_dec[9] <<  "\n";
     std::cout << a_dec[0] << a_dec[1] << a_dec[2] << a_dec[3] << a_dec[4] << a_dec[5] << a_dec[6] << a_dec[7] << a_dec[8] << a_dec[9] <<  "\n";
    std::cout << b_dec[0] << b_dec[1] << b_dec[2] << b_dec[3] << b_dec[4] << b_dec[5] << b_dec[6] << b_dec[7] << b_dec[8] << b_dec[9] <<  "\n";
   
    
/*
    // timing tests
    std::cout << "Tests\n";
    std::cout << time(0);
    
    std::cout << "End\n";
    

*/
    
    return 0;
}



