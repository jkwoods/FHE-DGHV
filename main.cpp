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



     // make pk
     //Pk pk_a = Pk::make_key(0); //toy
     
     //bool test = pk_a.assert_parameter_correctness();
     //std::cout << test << "\n";
     
    //std::vector<std::vector<int>> mata = { {1,0}, {0,1}};
    //std::vector<std::vector<int>> matb = { {1,0}, {1,0}};
    //Encoding aa = Encoding::matmul(mata, matb, pk_a);
    
    
     std::vector<int> z_vec = {0,0,0,0,0,0,1,0,0,0}; //,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

     Encoding zero = Encoding(pk_a, z_vec);


     std::vector<int> o_vec = {1,1,1,1,1,1,1,1,1,1}; //1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
     Encoding one = Encoding(pk_a, o_vec);
     
     //std::vector<int> a_vec = {1,1,1,1,1,0,0,0,0,0};
     //Encoding aa = Encoding::selector(a_vec,zero,one); //Encoding(pk_a, a_vec);
     
     //std::vector<int> b_vec = {0,1,0,1,0,1,0,1,0,1};
     //Encoding bb = zero.neg(); //Encoding(pk_a, b_vec);
     
     
     Encoding a = one+zero;
     Encoding b = one*zero;
     
     //std::vector<int> z_dec = aa.decode();

     std::vector<int> o_dec = b.decode();
     
     Encoding bad = b*zero;
     std::vector<int> a_dec = bad.decode();
     
     
     b.recode(1);
    

     Encoding good = b*one;
     std::vector<int> b_dec = good.decode();
     
     //print
     
     //std::cout << z_dec[0] << z_dec[1] << z_dec[2] << z_dec[3] << z_dec[4] << z_dec[5] << z_dec[6] << z_dec[7] << z_dec[8] << z_dec[9] <<  "\n";
     std::cout << o_dec[0] << o_dec[1] << o_dec[2] << o_dec[3] << o_dec[4] << o_dec[5] << o_dec[6] << o_dec[7] << o_dec[8] << o_dec[9] <<  "\n";
      std::cout << a_dec[0] << a_dec[1] << a_dec[2] << a_dec[3] << a_dec[4] << a_dec[5] << a_dec[6] << a_dec[7] << a_dec[8] << a_dec[9] <<  "\n";
     std::cout << b_dec[0] << b_dec[1] << b_dec[2] << b_dec[3] << b_dec[4] << b_dec[5] << b_dec[6] << b_dec[7] << b_dec[8] << b_dec[9] <<  "\n";
    
     
    
    /*
    
    //timing tests
    for(int i = 0; i < 1; i++){
        
        std::cout << "Key Size = " << i <<"\n";
        
        std::cout << "Making Keys\n";
        long start_time = time(0);
        Pk pk_a = Pk::make_key(i);
        long end_time = time(0);
        long run_time = end_time - start_time;
        std::cout << "key time: " << run_time <<  "\n";
        std::cout << "key correct? " << pk_a.assert_parameter_correctness() << "\n";
        
        int l = pk_a.p_l;
        std::vector<int> a_vec;
        std::vector<int> b_vec;
        srand(time(NULL));
        for(int k = 0; k < l; k++){
            int a_int = rand() % 2;
            a_vec.push_back(a_int);
            int b_int = rand() % 2;
            b_vec.push_back(b_int);
        }
        
        std::cout << "Encoding a\n";
        start_time = time(0);
        Encoding a = Encoding(pk_a, a_vec);
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "encode a time: " << run_time <<  "\n";
        
        std::cout << "Encoding b\n";
        start_time = time(0);
        Encoding b = Encoding(pk_a, b_vec);
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "encode b time: " << run_time <<  "\n";
        
        std::cout << "Adding a+b\n";
        start_time = time(0);
        Encoding c = a+b;
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "addition time: " << run_time <<  "\n";
        
        std::cout << "Mulitply a*b\n";
        start_time = time(0);
        Encoding d = a*b;
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "multiplication time: " << run_time <<  "\n";
        
        std::cout << "Recoding a*b\n";
        start_time = time(0);
        d.recode();
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "recode time: " << run_time <<  "\n";
        
        std::cout << "Decoding a\n";
        start_time = time(0);
        std::vector<int> a_dec = a.decode();
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "decode a time: " << run_time <<  "\n";
        
        std::cout << "Decoding c\n";
        start_time = time(0);
        std::vector<int> c_dec = c.decode();
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "decode c time: " << run_time <<  "\n";
        
        std::cout << "Decoding d\n";
        start_time = time(0);
        std::vector<int> d_dec = d.decode();
        end_time = time(0);
        run_time = end_time - start_time;
        std::cout << "decode d time: " << run_time <<  "\n";
        
        std::cout << "Tests Done\n";
        
        
        
        
        
        
    }
    
    
    */
    //bool test = pk_a.assert_parameter_correctness();
    
    
    
    
    
    //Encoding one = Encoding(pk_a, o_vec);
    
    //std::vector<int> a_vec = {1,1,1,1,1,0,0,0,0,0};
    //Encoding aa = Encoding::selector(a_vec,zero,one); //Encoding(pk_a, a_vec);
    
    //std::vector<int> b_vec = {0,1,0,1,0,1,0,1,0,1};
    //Encoding bb = zero.neg(); //Encoding(pk_a, b_vec);
    
    
    //Encoding a = one+zero;
    //Encoding b = one*zero;
    
    //std::vector<int> z_dec = a.decode();
    //std::vector<int> o_dec = b.decode();
    
    //Encoding bad = b*zero;
    //std::vector<int> a_dec = bad.decode();
    
    /*
    b.recode();
    Encoding good = b*zero;
    std::vector<int> b_dec = good.decode();
    
    //print
    
    std::cout << z_dec[0] << z_dec[1] << z_dec[2] << z_dec[3] << z_dec[4] << z_dec[5] << z_dec[6] << z_dec[7] << z_dec[8] << z_dec[9] <<  "\n";
    std::cout << o_dec[0] << o_dec[1] << o_dec[2] << o_dec[3] << o_dec[4] << o_dec[5] << o_dec[6] << o_dec[7] << o_dec[8] << o_dec[9] <<  "\n";
     std::cout << a_dec[0] << a_dec[1] << a_dec[2] << a_dec[3] << a_dec[4] << a_dec[5] << a_dec[6] << a_dec[7] << a_dec[8] << a_dec[9] <<  "\n";
    std::cout << b_dec[0] << b_dec[1] << b_dec[2] << b_dec[3] << b_dec[4] << b_dec[5] << b_dec[6] << b_dec[7] << b_dec[8] << b_dec[9] <<  "\n";
   
     */
    
    return 0;
}



