//
//  Deltas.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "Deltas.hpp"
#include "Pk.hpp"
#include <cmath>

Deltas::Deltas(Pk pk, int lenv, int rho, int cr)
    : r_pk(pk), r_lenv(lenv), r_rho(rho), r_cr(cr), r_pri(makePri()), r_deltas(), r_Chi(r_pri.r_list) {
        
        
        makeDeltas();
        makeDeltaList();
        makeState();
    }

Deltas::~Deltas(){ //TODO -vectors??
    for(int i = 0; i < r_x.size(); i++) { mpz_clear(r_x[i]); }
    for(int i = 0; i < r_deltas.size(); i++) { mpz_clear(r_deltas[i]); }
    for(int i = 0; i < r_Chi.size(); i++) { mpz_clear(r_Chi[i]); }

}

void Deltas::makeDeltaList(){
    
    for(int i = 0; i < r_pri.r_len; i++){
        mpz_init(r_x[i]);
        mpz_sub(r_x[i], r_Chi[i], r_deltas[i]);
        
        //x.push_back(r_Chi[i]-r_deltas[i]);
    }
}

PseudoRandomInts Deltas::makePri(){
    PseudoRandomInts pri = PseudoRandomInts(r_pk.p_x0, r_lenv);
    return pri;
}

void Deltas::makeDeltas(){
    //make deltas
    std::vector<std::vector<mpz_t>> r;
    std::vector<mpz_t> E;
    
    mpz_t e_help;
    //int e_help = pow(2,(r_pk.p_lam+r_pk.p_logl+(r_pk.p_l*r_pk.p_eta)));
    
    mpz_init_set_ui(e_help, r_pk.p_eta);
    mpz_mul_ui(e_help, e_help, r_pk.p_l);
    mpz_add_ui(e_help, e_help, r_pk.p_logl);
    mpz_add_ui(e_help, e_help, r_pk.p_lam);
    
    for(int i = 0; i < r_lenv; i++){
        for(int j = 0; j < r_pk.p_l; j++){
            mpz_init(r[i][j]);
            
            random_element_pow2(r[i][j], r_rho+1, r_rho, d_state); //-2 and 2
            //int rand = random_element(pow(-2,r_rho+1), pow(2,r_rho));
        }
        mpz_init(E[i]);

        mpz_urandomb(E[i], state, e_help); //2^e_help
        
        //E.push_back(Ei / r_pk.p_pi);
        mpz_fdiv_q(E[i], E[i], r_pk.p_pi);
    }
    
    std::vector<mpz_t> crts;
    if (r_cr==0){ //x
        for(int i = 0; i < r_lenv; i++){
            mpz_init(crts[i]);
            std::vector<mpz_t> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                mpz_init(crt_term[i]);
                mpz_mul_ui(crt_term[i], r[i][j], 2);
                
                //crt_term.push_back(2*r[i][j]);
            }
            CRT(crts[i], r_pk.p_p, crt_term);
            for (int j = 0; j < r_pk.p_l; j++){ mpz_clear(crt_term[j]); }
            //crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else if (r_cr==1){ //xi
        for(int i = 0; i < r_lenv; i++){
            mpz_init(crts[i]);
            std::vector<mpz_t> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                mpz_init(crt_term[i]);
                mpz_mul_ui(crt_term[i], r[i][j], 2);
                mpz_add_ui(crt_term[i], crt_term[i], kd(j,i));
                
                //crt_term.push_back(2*r[i][j]+kd(j,i));
            }
            CRT(crts[i], r_pk.p_p, crt_term);
            for (int j = 0; j < r_pk.p_l; j++){ mpz_clear(crt_term[j]); }
            //crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else if (r_cr==2){ //ii
        for(int i = 0; i < r_lenv; i++){
            mpz_init(crts[i]);
            std::vector<mpz_t> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                mpz_init(crt_term[i]);
                mpz_mul_ui(crt_term[i], r[i][j], 2);
                mpz_add_ui(crt_term[i], crt_term[i], (kd(j,i)*(pow(2,r_pk.p_rhoi+1))));
                
                //crt_term.push_back(2*r[i][j]+(kd(j,i)*(pow(2,r_pk.p_rhoi+1))));
            }
            CRT(crts[i], r_pk.p_p, crt_term);
            for (int j = 0; j < r_pk.p_l; j++){ mpz_clear(crt_term[j]); }
            //crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    } else { //o
        for(int i = 0; i < r_lenv; i++){
            mpz_init(crts[i]);
            std::vector<mpz_t> crt_term;
            for (int j = 0; j < r_pk.p_l; j++){
                mpz_init(crt_term[i]);
                mpz_mul_ui(crt_term[i], r[i][j], 2);
                mpz_add_ui(crt_term[i], crt_term[i], r_pk.p_vert_s[i][j]);
                
                //crt_term.push_back(2*r[i][j]+r_pk.p_vert_s[i][j]);
            }
            CRT(crts[i], r_pk.p_p, crt_term);
            for (int j = 0; j < r_pk.p_l; j++){ mpz_clear(crt_term[j]); }
            //crts.push_back(CRT(r_pk.p_p, crt_term));
        }
    }
    
    for (int i = 0; i < r_lenv; i++){
        mpz_init(r_deltas[i]);
        
        mpz_t chi_temp;
        mpz_init(chi_temp);
        mpz_mod(chi_temp, r_Chi[i], r_pk.p_pi);
        //int chi_temp = mod(r_Chi[i],r_pk.p_pi);
        
        mpz_mul(r_deltas[i], E[i], r_pk.p_pi);
        mpz_add(r_deltas[i], r_deltas[i], chi_temp);
        mpz_sub(r_deltas[i], r_deltas[i], crts[i]);
        //int delta_term = chi_temp+(E[i]*r_pk.p_pi)-crts[i];
        
    }

    for(int i = 0; i < r_lenv; i++){
        for(int j = 0; j < r_pk.p_l; j++){ mpz_clear(r[i][j]); }
        mpz_clear(E[i]);
    }
    mpz_clear(e_help);
}

void Deltas::makeState(){
    gmp_randinit_mt(d_state);
    gmp_randseed_ui(d_state, time(0)); //time as seed TODO - check this shit
}
