//
//  Pk.hpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#ifndef Pk_hpp
#define Pk_hpp

#include <stdio.h>
#include <vector>
#include "utils.hpp"

class Pk{
private:
    int p_lam;
    int p_rho;
    int p_rhoi;
    int p_eta;
    int p_gam;
    int p_Theta;
    int p_theta;
    int p_n;
    int p_kap;
    int p_alpha;
    int p_alphai;
    int p_tau;
    int p_l;
    int p_logl;
    std::vector<int> p_p;
    int p_pi;
    int p_q0;
    int p_x0;
    std::vector<int> p_x;
    std::vector<int> p_xi;
    std::vector<int> p_ii;
    int p_B;
    std::vector<std::vector<int>> p_s;
    std::vector<std::vector<int>> p_vert_s;
    std::vector<int> p_u;
    std::vector<int> p_o;
    
    //helper
    std::vector<int> make_p();
    int make_pi();
    int make_q0();
    std::vector<int> make_x();
    std::vector<int> make_xi();
    std::vector<int> make_ii();
    std::vector<std::vector<int>> make_s();
    std::vector<std::vector<int>> make_vert_s();
    std::vector<int> make_u();
    std::vector<int> make_o();
    
public:
    Pk(int lam, int rho, int rhoi, int eta, int gam, int Theta, int theta, int kap, int alpha, int alphai, int tau, int l, int n=4);
    int encode(std::vector<int> m);
    std::vector<int> decode(int c);
    std::vector<int> decode_squashed(int c);
    int recode(int c);
    int H_add(int c1, int c2);
    int H_mult(int c1, int c2);

};


#endif /* Pk_hpp */
