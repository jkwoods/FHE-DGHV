//
//  Deltas.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "Deltas.hpp"


Deltas::Deltas(Pk pk, int lenv, int rho, int cr)
    : r_pk(pk), r_lenv(lenv), r_rho(rho), r_cr(cr), r_pri(PseudoRandomInts(pk.p_x0, lenv)), r_deltas(lenv), r_Chi(lenv), r_x(lenv) {
        
        makeDeltas();
        makeDeltaList(); //'x'
        
        makeState(); //TODO?
    }

Deltas::~Deltas(){
    //TODO -vectors?
}
    
void Deltas::makeDeltaList(){
    for(int i = 0; i < r_pri.r_len; i++){
        r_x[i] = r_Chi[i]-r_deltas[i];
    }
}

void Deltas::makeDeltas(){
    //make deltas
    std::vector<std::vector<mpz_class>> r(r_lenv, std::vector<mpz_class> (r_pk.p_l)); //correct dim?
    std::vector<mpz_class> E(r_lenv);
    
    mpz_class e_help = pow(2,(r_pk.p_lam+r_pk.p_logl+(r_pk.p_l*r_pk.p_eta)));
    
    for(int i = 0; i < r_lenv; i++){
        for(int j = 0; j < r_pk.p_l; j++){
            r[i][j] = random_element(pow(-2,r_rho+1), pow(2,r_rho));
        }
        E[i] = e_help / r_pk.p_pi; //floor
    }
    
    std::vector<mpz_class> crts(r_lenv);
    if (r_cr==0){ //x
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[i] = (2*r[i][j]);
            }
            crts[i] = CRT(r_pk.p_p, crt_term);

        }
    } else if (r_cr==1){ //xi
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[i] = (2*r[i][j]+kd(j,i));
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    } else if (r_cr==2){ //ii
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[i] = (2*r[i][j]+(kd(j,i)*(pow(2,r_pk.p_rhoi+1))));
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    } else { //o
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[i] = (2*r[i][j]+r_pk.p_vert_s[i][j]);
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    }
    
    for (int i = 0; i < r_lenv; i++){
        mpz_class chi_temp = mod(r_Chi[i],r_pk.p_pi);
        
        r_deltas[i] = chi_temp+(E[i]*r_pk.p_pi)-crts[i];
    }
}

void Deltas::makeState(){
    gmp_randinit_mt(d_state);
    gmp_randseed_ui(d_state, time(0)); //time as seed TODO - check this shit
}
