//
//  Deltas.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  Copyright © 2019 Woods, Jess. All rights reserved.
//

#include "Deltas.hpp"
#include "Pk.hpp"
#include <cmath>

Deltas::Deltas(Pk pk, int lenv, int rho, int cr)
    : r_pk(pk), r_lenv(lenv), r_rho(rho), r_cr(cr), r_pri(makePri()), r_deltas(makeDeltas()), r_Chi(r_pri.getList()) {}

std::vector<int> Deltas::getDeltaList(){
    std::vector<int> x;
    for(int i = 0; i < r_pri.r_len; i++){
        x.push_back(r_Chi[i]-r_deltas[i]);
    }
    return x;
}

PseudoRandomInts Deltas::makePri(){
    PseudoRandomInts pri = PseudoRandomInts(r_pk.p_x0, r_lenv);
    return pri;
}

std::vector<int> Deltas::makeDeltas(){
    //make deltas
    std::vector<std::vector<int>> r;
    std::vector<int> E;
    for(int i = 0; i < r_lenv; i++){
        std::vector<int> ri;
        for(int i = 0; i < r_pk.p_l; i++){
            int rand = random_element(pow(-2,r_rho+1), pow(2,r_rho));
            ri.push_back(rand);
        }
        r.push_back(ri);
        
        int u = pow(2,(r_pk.p_lam+r_pk.p_logl+(r_pk.p_l*r_pk.p_eta)));
        int Ei = random_element(0, u);
        E.push_back(Ei / r_pk.p_pi);
    }
    
    std::vector<int> crts;
    if (r_cr==0){ //x
        for(int i = 0; i < r_lenv; i++){
            std::vector<int> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term.push_back(2*r[i][j]);
            }
            crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else if (r_cr==1){ //xi
        for(int i = 0; i < r_lenv; i++){
            std::vector<int> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term.push_back(2*r[i][j]+kd(j,i));
            }
            crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else if (r_cr==2){ //ii
        for(int i = 0; i < r_lenv; i++){
            std::vector<int> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term.push_back(2*r[i][j]+(kd(j,i)*(pow(2,r_pk.p_rhoi+1))));
            }
            crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else { //o
        for(int i = 0; i < r_lenv; i++){
            std::vector<int> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term.push_back(2*r[i][j]+r_pk.p_vert_s[i][j]);
            }
            crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    }
    
    std::vector<int> deltas;
    for (int i = 0; i < r_lenv; i++){
        int chi_temp = mod(r_Chi[i],r_pk.p_pi);
        int delta_term = chi_temp+(E[i]*r_pk.p_pi)-crts[i];
        deltas.push_back(delta_term);
    }

    return deltas;
}
