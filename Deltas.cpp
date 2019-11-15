//
//  Deltas.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "Deltas.hpp"


Deltas::Deltas(Pk& pk, int lenv, int rho, int cr)
: r_pk(pk), r_rho(rho), r_cr(cr), r_deltas(lenv), r_Chi(lenv), r_lenv(lenv), r_pri(PseudoRandomInts(r_pk.p_x0, lenv)), r_x(lenv)
{
        
        makeDeltas();
        makeDeltaList(); //'x'
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
    //make class state
    gmp_randclass p_class_state (gmp_randinit_mt);
    p_class_state.seed(time(0)); //TODO
    
    
    //make deltas
    std::vector<std::vector<mpz_class>> r(r_lenv, std::vector<mpz_class> (r_pk.p_l)); //correct dim?
    std::vector<mpz_class> E(r_lenv);
    
    mpz_class e_help = power(2,(r_pk.p_lam+r_pk.p_logl+(r_pk.p_l*r_pk.p_eta))); //is this contained by int? TODO
    
    for(int i = 0; i < r_lenv; i++){
        for(int j = 0; j < r_pk.p_l; j++){
            mpz_class lb = power(-2,r_rho+1);
            mpz_class ub = power(2,r_rho);
            r[i][j] = p_class_state.get_z_range(ub-lb);
            r[i][j] = r[i][j] + lb;
        }
        mpz_class rand = p_class_state.get_z_range(e_help);
        E[i] = rand / r_pk.p_pi; //floor
    }
    
    std::vector<mpz_class> crts(r_lenv);
    if (r_cr==0){ //x
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[j] = 2*r[i][j];
            }
            crts[i] = CRT(r_pk.p_p, crt_term);

        }
    } else if (r_cr==1){ //xi
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[j] = (2*r[i][j]+kd(j,i));
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    } else if (r_cr==2){ //ii
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[j] = (2*r[i][j]+(kd(j,i)*(power(2,r_pk.p_rhoi+1))));
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    } else { //o
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(r_pk.p_l);
            for (int j = 0; j < r_pk.p_l; j++){
                crt_term[j] = (2*r[i][j]+r_pk.p_vert_s[i][j]);
            }
            crts[i] = CRT(r_pk.p_p, crt_term);
        }
    }
    
    for (int i = 0; i < r_lenv; i++){
        mpz_class chi_temp = r_Chi[i] % r_pk.p_pi;
        
        r_deltas[i] = chi_temp+(E[i] * r_pk.p_pi)-crts[i];
    }
}
