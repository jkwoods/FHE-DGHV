
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
#include <chrono>

#ifdef _OPENMP
    #include "omp.h"
#else
    #define omp_get_wtime() std::clock()
#endif

typedef std::chrono::high_resolution_clock Clock;

/* TODO
- Optimize recode
 

 */

int main(int argc, const char * argv[]) {


     std::cout << "Making Keys\n";
    // Pk pk_a = Pk(lam, rho, eta, gam, Theta, alpha, tau, l);

//
     // make pk
    Pk pk_a = Pk::make_key(0); 
    double start = omp_get_wtime();
    //auto s = Clock::now();
    pk_a = Pk::make_key(0); //toy
    double end = omp_get_wtime();
    //auto e = Clock::now();       
    //std::cout << "key time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(e - s).count()  << "\n";
    std::cout << "key time: " << (end - start) << "\n";
 
     bool test = pk_a.assert_parameter_correctness();
     std::cout << test << "\n";
     
    //std::vector<std::vector<int>> mata = { {1,0}, {0,1}};
    //std::vector<std::vector<int>> matb = { {1,0}, {1,0}};
    //Encoding aa = Encoding::matmul(mata, matb, pk_a);

    //bool test = pk_a.assert_parameter_correctness();
    //std::cout << test << "\n";
    
    
    std::vector<int> z_vec = {1,0,1,0,1,0,1,0,1,0}; //,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    Encoding zero = Encoding(pk_a, z_vec);
    start = omp_get_wtime();
    //s = Clock::now();
    zero = Encoding(pk_a, z_vec);
    //e = Clock::now();
    end = omp_get_wtime();    
    //std::cout << "encoding time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(e - s).count() << "\n";
    std::cout << "encoding time: " << (end - start) << "\n";

    std::vector<int> o_vec = {1,1,1,1,1,1,1,1,1,1}; //1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
    Encoding one = Encoding(pk_a, o_vec);
    
    std::vector<int> a_vec = {1,1,1,1,1,0,0,0,0,0};
    Encoding aa = Encoding::selector(a_vec,zero,one); //Encoding(pk_a, a_vec);
    
    std::vector<int> b_vec = {0,1,0,1,0,1,0,1,0,1};
    Encoding bb = zero.neg(); //Encoding(pk_a, b_vec);
    
    //start = omp_get_wtime();
    Encoding a = one+zero;
    //end = omp_get_wtime();    
    //std::cout << "addition time: " << (end-start) << "\n";

    //start = omp_get_wtime();
    Encoding b = one*zero;
    //end = omp_get_wtime();
    //std::cout << "mult time: " << (end-start) << "\n";    

    std::vector<int> z_dec = zero.decode();
    start = omp_get_wtime();
    //s = Clock::now();
    z_dec = zero.decode();
    //e = Clock::now();
    end = omp_get_wtime();
    //std::cout << "decoding time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(e - s).count() << "\n";    
    std::cout << "decoding time: " << (end - start) << "\n";


    //start = omp_get_wtime();
    std::vector<int> o_dec = one.decode();
    end = omp_get_wtime();
    //std::cout << "decoding time ones: " << (end-start) << "\n";
    
    //start = omp_get_wtime();
    std::vector<int> a_dec = (Encoding(pk_a, 1)).decode();
    //end = omp_get_wtime();
    //std::cout << "decoding with noise: " << (end-start) << "\n";

    
    //start = omp_get_wtime();
    a.recode(0);
    //end = omp_get_wtime();
    //std::cout << "recoding time: " << (end-start) << "\n";

    start = omp_get_wtime();
    //s = Clock::now();
    b.recode(0);
    //e = Clock::now();
    end = omp_get_wtime();
    //std::cout << "recoding time: " << std::chrono::duration_cast<std::chrono::nanoseconds>(e - s).count() << "\n";
    std::cout << "recoding time: " << (end - start) << "\n";

    //start = omp_get_wtime();
    //std::vector<int> b_dec = b.decode();
    //end = omp_get_wtime();
    //std::cout << "decoding after recode: " << (end-start) << "\n";

    Encoding good = b*zero;
    std::vector<int> b_dec = b.decode();
    
    //print
    
    std::cout << z_dec[0] << z_dec[1] << z_dec[2] << z_dec[3] << z_dec[4] << z_dec[5] << z_dec[6] << z_dec[7] << z_dec[8] << z_dec[9] <<  "\n";
    std::cout << o_dec[0] << o_dec[1] << o_dec[2] << o_dec[3] << o_dec[4] << o_dec[5] << o_dec[6] << o_dec[7] << o_dec[8] << o_dec[9] <<  "\n";
     std::cout << a_dec[0] << a_dec[1] << a_dec[2] << a_dec[3] << a_dec[4] << a_dec[5] << a_dec[6] << a_dec[7] << a_dec[8] << a_dec[9] <<  "\n";
    std::cout << b_dec[0] << b_dec[1] << b_dec[2] << b_dec[3] << b_dec[4] << b_dec[5] << b_dec[6] << b_dec[7] << b_dec[8] << b_dec[9] <<  "\n";
   
    return 0;
}



