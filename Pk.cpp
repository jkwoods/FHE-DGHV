//
//  Pk.cpp
//  FHE
//
//  Created by Woods, Jess on 10/15/19.
//  work done at Oak Ridge National Lab
//


// TODO
// put the y/u/z etc stuff in the correct places
// pre-generate a bunch of random ## for encoding process
// make longs??


#include "Pk.hpp"
#include <iostream>
#include "utils.hpp"

Pk::Pk(int lam, int rho, int eta, int gam, int Theta, int alpha, int tau, int l, int n)
: p_lam(lam), p_rho(rho), p_rhoi(rho+lam), p_eta(eta), p_gam(gam), p_Theta(Theta), p_theta(Theta/l), p_kap(64*(gam/64+1)-1), p_alpha(alpha), p_alphai(alpha+lam), p_tau(tau), p_l(l), p_logl(int (round(log2(l)))), p_p(l), p_pi(1), p_q0(1), p_x0(1), p_seed_x(time(0)), p_seed_xi(time(0)+3), p_seed_ii(time(0)+8), p_seed_o(time(0)+19), p_seed_op(time(0)+27), p_seed_u(time(0)+33), p_x_D(tau), p_xi_D(l), p_ii_D(l), p_n(n), p_B(l), p_s(l,std::vector<int>(Theta)), p_vert_s(Theta,std::vector<int>(l)), p_u_front(p_l), p_o_D(Theta), p_op_D(Theta) //for kap, c++ trucates (rounds down) so its all gud
{
    
    std::cout << "Parameters secure and correct? " << this->assert_parameter_correctness() << "\n";
    
    //make t state
    gmp_randstate_t p_t_state;
    gmp_randinit_mt(p_t_state);
    gmp_randseed_ui(p_t_state, time(0)); //time is not great - better TODO
    
    make_p(p_t_state);
    make_pi();
    make_q0(p_t_state);
    make_x0();
    make_x_Delta();
    make_xi_Delta();
    make_ii_Delta();
    make_s();
    make_vert_s();
    make_u();
    make_y();
    make_o_Delta();
    make_op_Delta();
    
}

//Pk::~Pk(){
    // TODO
    
    //state clear
//    gmp_randclear(p_t_state);
//}

mpz_class Pk::encode(std::vector<int> m){
    //recover x, xi, ii
    PseudoRandomInts chi_x = PseudoRandomInts(p_x0, p_tau, p_seed_x);
    PseudoRandomInts chi_xi = PseudoRandomInts(p_x0, p_l, p_seed_xi);
    PseudoRandomInts chi_ii = PseudoRandomInts(p_x0, p_l, p_seed_ii);
    
    std::vector<mpz_class> x = make_x_list(chi_x, p_x_D);
    std::vector<mpz_class> xi = make_x_list(chi_xi, p_xi_D);
    std::vector<mpz_class> ii = make_x_list(chi_ii, p_ii_D);
    
    
    //make class state
    gmp_randclass p_class_state (gmp_randinit_mt);
    p_class_state.seed(time(0)); //TODO
    
    //m*xi
    std::vector<mpz_class> m_xi(p_l);
    for (int i = 0; i < p_l; i++){
        m_xi[i] = m[i]*xi[i];
    }
    
    //bi*ii
   
    std::vector<mpz_class> bi_ii(p_l);
    for (int i = 0; i < p_l; i++){
        
        mpz_class lb = power(-2,p_alphai);
        mpz_class ub = power(2,p_alphai);
        mpz_class bi = p_class_state.get_z_range(ub-lb);
        bi = bi + lb;
        bi_ii[i] = bi*ii[i];
    }
    
    //b*x
    std::vector<mpz_class> b_x(p_tau);
    for (int i = 0; i < p_tau; i++){

        mpz_class lb = power(-2,p_alpha);
        mpz_class ub = power(2,p_alpha);
        mpz_class b = p_class_state.get_z_range(ub-lb);
        b = b + lb;
        
        b_x[i] = b*x[i];
    }
    
    //summation
    mpz_class bigsum = sum_array(m_xi) + sum_array(bi_ii) + sum_array(b_x);
    
    mpz_class c = modNear(bigsum, p_x0);
    
    return c;

}

std::vector<int> Pk::decode(mpz_class c){
    

    std::vector<int> m(p_l);
    for (int i = 0; i < p_l; i++){
        //std::cout << "slot " << i << "\n";
        //std::cout << "pi= " << p_p[i] << "\n";
        
        mpz_class mn = modNear(c,p_p[i]);
        //std::cout << "modNear " << mn << "\n";
        
        mpz_class conv = floor_mod(mn,2);
        //std::cout << "mod2 " << conv << "\n";
        
        int i_conv = (int) conv.get_si(); //hopefully right
        m[i] = (i_conv);
    }
    return m;

}

std::vector<int> Pk::decode_squashed(mpz_class c){ //TODO gen
    std::vector<int> temp;
    return temp;
}

std::vector<std::vector<int>> Pk::expand(mpz_class c){
    //recover u
    
    mpz_class pwr = power(2,p_kap+1);
    PseudoRandomInts priu = PseudoRandomInts(pwr, p_Theta, p_seed_u);
             
    std::vector<mpz_class> u = make_full_u(priu);
    
    //std::vector<mpq_class> y = make_y();
    
    std::vector<std::vector<int>> z(p_Theta,std::vector<int>(p_n+1));
    
    for (int i = 0; i < p_Theta; i++){
        /*
        mpq_class zi = c * y[i];
        mpq_class mod = mod_2_f(zi);
        //std::cout << mpf_class(mod) << "\n";
        
        //convert each zi to vector of binary
        //mult by 2^n (bits of precision)
        mpq_class zi_mult = mod * (pow(2,p_n));
        //std::cout << mpf_class(zi_mult) << "\n";
        
        //convert to int (cut off ends)
        mpz_class zi_conv = mpz_class(zi_mult);
        if ((zi_mult - zi_conv) > 0.5){
            zi_conv += 1;
        }
        //std::cout << zi_conv << "\n";
        
        if (zi_conv.fits_sint_p() == 0){ //doesn't fit
            std::cout << "Error in generation of z\n";
        }
        int zi_round = zi_conv.get_si(); //round
        //std::cout << zi_round << "\n";
        
        //put in binary
        std::vector<int> ex_z(p_n+1);
        ex_z = to_binary(zi_round, p_n+1); //[lsb, ...., msb] max 31 (can we ever possibly get a 32? what then) TODO
        //print_vec(ex_z);
        
        z[i] = ex_z;
         */
        
       
    
        unsigned long zi = directScal(p_kap, c.get_mpz_t(), u[i].get_mpz_t());
        int new_zi = zi >> ((32-(p_n+2))+1) >> 1;
        
        //put in binary
        std::vector<int> ex_z(p_n+1);
        ex_z = to_binary(new_zi, p_n+1); //[lsb, ...., msb] max 31
        //print_vec(ex_z);
        
        z[i] = ex_z;
    }
    return z;
}

mpz_class Pk::recode_work(mpz_class c, std::vector<mpz_class> o){
    std::vector<std::vector<int>> z = expand(c);
    
    std::vector<std::vector<mpz_class>> z_mult(p_Theta,std::vector<mpz_class>(p_n+1));
    
    //z * sk - correct
    for(int i = 0; i < p_Theta; i++){
        for(int j = 0; j < p_n+1; j++){
            z_mult[i][j] = (z[i][j] * o[i]);
        }
    }
    
    //add to make as
    std::vector<mpz_class> a(p_n+1);
    for(int i = 0; i < p_Theta; i++){
        a = sum_binary(a,z_mult[i]);

        //for (int j = 0; j < p_n+1; j++){
            //a[j] = floor_mod(a[j], p_x0);
            //std::cout << "a; i= " << i << "\n";
            //print_vec(decode(a[j]));
        //}
        
    }
    
    //print_vec(decode(a[a.size()-2]));
    //std::cout << "c & 1= " << (c & 1) << "\n";
    
    mpz_class two = a[a.size()-1] + a[a.size()-2]; //correct??
    mpz_class c_prime = two + (c & 1);
    
    return c_prime;
}

mpz_class Pk::recode(mpz_class c){
    //recover o
    PseudoRandomInts chi_o = PseudoRandomInts(p_x0, p_Theta, p_seed_o);
    std::vector<mpz_class> o = make_x_list(chi_o, p_o_D);
    
    return recode_work(c, o);
}

mpz_class Pk::recode_and_permute(mpz_class c){
    //recover permuted o
    PseudoRandomInts chi_op = PseudoRandomInts(p_x0, p_Theta, p_seed_op);
    std::vector<mpz_class> op = make_x_list(chi_op, p_op_D);
   
    return recode_work(c, op);
}

mpz_class Pk::H_add(mpz_class c1, mpz_class c2){
    mpz_class c = floor_mod((c1+c2),p_x0);
    return c;
}

mpz_class Pk::H_mult(mpz_class c1, mpz_class c2){
    mpz_class c = floor_mod((c1*c2),p_x0);
    return c;
}


//PRIVATE HELPER
void Pk::make_p(gmp_randstate_t p_t_state){
    for (int i = 0; i < p_l; i++){
       p_p[i] = random_prime_w(p_eta, p_t_state); //weird range 2^(n-1), 2^n
    }
}

void Pk::make_pi(){ //prod of all p[i]
    p_pi = 1;
    for (int i = 0; i < p_l; i++){
        p_pi = p_pi*p_p[i];
    }
}

void Pk::make_q0(gmp_randstate_t p_t_state){
    p_q0 = power(2,p_gam);
    mpz_class comp = floor_div(p_q0,p_pi);
    
    while (p_q0 > comp){  
        mpz_class q0_prime1 = random_prime_f0(pow(p_lam,2), p_t_state); //2^entry //need to be long?TODO
        mpz_class q0_prime2 = random_prime_f0(pow(p_lam,2), p_t_state); //2^entry
        
        p_q0 = q0_prime1*q0_prime2;
    }
}

void Pk::make_x0(){
    p_x0 = p_pi * p_q0;
}

void Pk::make_x_Delta(){ //TODO DELTAS - TODO initialize list
    PseudoRandomInts pri = PseudoRandomInts(p_x0, p_tau, p_seed_x);
    p_x_D = make_Deltas(pri, p_tau, p_rhoi-1, 0);
}

void Pk::make_xi_Delta(){
    PseudoRandomInts pri = PseudoRandomInts(p_x0, p_l, p_seed_xi);
    p_xi_D = make_Deltas(pri, p_l, p_rho, 1);
}

void Pk::make_ii_Delta(){
    PseudoRandomInts pri = PseudoRandomInts(p_x0, p_l, p_seed_ii);
    p_ii_D = make_Deltas(pri, p_l, p_rho, 2);
}

std::vector<mpz_class> Pk::make_x_list(PseudoRandomInts chi, std::vector<mpz_class> deltas){
    std::vector<mpz_class> x(chi.r_len);
    for(int i = 0; i < chi.r_len; i++){
        x[i] = chi.r_list[i]-deltas[i];
    }
    return x;
}

void Pk::make_s(){
    for(int j = 0; j < p_l; j++){
        //all initialized to 0 originally
            p_s[j][j] = 1;
    }
    
    for(int t = 1; t < p_theta; t++){
        std::vector<int> sri = random_sample(p_B, p_l);
        for(int j = 0; j < p_l; j++){
            int k = (p_B*t)+sri[j];
            p_s[j][k] = 1;
        }
    }
    
    //print
    /*
    for (int j = 0; j < 15; j++)
    {
        for (int i = 0; i < p_l; i++)
        {
            std::cout << p_s[i][(j*10)+0] << p_s[i][(j*10)+1] << p_s[i][(j*10)+2] << p_s[i][(j*10)+3] << p_s[i][(j*10)+4] << p_s[i][(j*10)+5] << p_s[i][(j*10)+6] << p_s[i][(j*10)+7] << p_s[i][(j*10)+8] << p_s[i][(j*10)+9] << "\n";
        }
        std::cout << "\n\n";
    }
    */
}

void Pk::make_vert_s(){
    for(int i = 0; i < p_Theta; i++){
        for(int j = 0; j < p_l; j++){
            p_vert_s[i][j] = p_s[j][i];
        }
    }
}


void Pk::make_u(){ //in making Pk
    mpz_class pwr = power(2,p_kap+1);
    PseudoRandomInts priu = PseudoRandomInts(pwr, p_Theta, p_seed_u);
    p_u_front = make_front_U(priu);
}

std::vector<mpq_class> Pk::make_y(){ //in recode
    std::vector<mpq_class> y(p_Theta); //rational
    
    mpz_class pwr = power(2,p_kap+1);
    PseudoRandomInts priu = PseudoRandomInts(pwr, p_Theta, p_seed_u);
    
    std::vector<mpz_class> u = make_full_u(priu);
    
    mpz_class div = power(2,p_kap);
       for (int i = 0; i < p_Theta; i++){
           y[i] = mpq_class(u[i],div); //rational u[i]/(2^kap)
           y[i].canonicalize();
       }
    return y;
}

std::vector<mpz_class> Pk::make_full_u(PseudoRandomInts priu){
    std::vector<mpz_class> full_u(p_Theta);
    
    for (int i = 0; i < p_Theta; i++){
        if (i < p_l){
            full_u[i] = p_u_front[i];
        } else {
            full_u[i] = priu.r_list[i];
        }
    }
    return full_u;
}

void Pk::make_o_Delta(){
    PseudoRandomInts pri = PseudoRandomInts(p_x0, p_Theta, p_seed_o);
    p_o_D = make_Deltas(pri, p_Theta, p_rho, 3);
}

void Pk::make_op_Delta(){
    PseudoRandomInts pri = PseudoRandomInts(p_x0, p_Theta, p_seed_op);
    p_op_D = make_Deltas(pri, p_Theta, p_rho, 4);
}

bool Pk::assert_parameter_correctness(){
    bool a = p_rho >= 2*p_lam; //brute force noise attack
    bool b = p_eta >= p_alphai + p_rhoi + 1 + log2(p_l); // correct decoding
    bool c = p_eta >= p_rho * (p_lam*(pow(log(p_lam),2))); //squashed decode circut
    //bool d = p_gam > pow(p_eta, 2) * log(p_lam); //lattice attack
    
    bool e = (p_alpha * p_tau) >= p_gam + p_lam; //leftover hash lemma
    bool f = p_tau >= p_l * (p_rhoi + 2) + p_lam; //leftover hash lemma
    bool g = (p_Theta % p_l == 0);
    
    return a && b && c && e && f && g;
}

Pk Pk::make_key(int size){
    int lam = 12;
    int rho = 26;
    int eta = 15256;
    int gam = 147456;
    int Theta = 150;
    int alpha = 936;
    int tau = 742;
    int l = 10;
    
    return Pk(lam, rho, eta, gam, Theta, alpha, tau, l);
}

//old deltas
std::vector<mpz_class> Pk::make_Deltas(PseudoRandomInts r_chi, int r_lenv, int r_rho, int r_cr){
    std::cout << "r_lenv=" << r_lenv << "\n";
    
    std::vector<mpz_class> deltas(r_lenv);
    
    //make class state
    gmp_randclass p_class_state (gmp_randinit_mt);
    p_class_state.seed(time(0)); //TODO
    
    
    //make deltas
    std::vector<std::vector<mpz_class>> r(r_lenv, std::vector<mpz_class> (p_l)); //correct dim?
    std::vector<mpz_class> E(r_lenv);
    
    mpz_class e_help = power(2,(p_lam+p_logl+(p_l*p_eta))); //is this contained by int? TODO
    
    for(int i = 0; i < r_lenv; i++){
        for(int j = 0; j < p_l; j++){
            mpz_class lb = power(-2,r_rho+1);
            mpz_class ub = power(2,r_rho);
            r[i][j] = p_class_state.get_z_range(ub-lb);
            r[i][j] = r[i][j] + lb;
        }
        mpz_class rand = p_class_state.get_z_range(e_help);
        E[i] = floor_div(rand,p_pi); //floor
    }
    
    std::vector<mpz_class> crts(r_lenv);
    if (r_cr==0){ //x
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(p_l);
            for (int j = 0; j < p_l; j++){
                crt_term[j] = 2*r[i][j];
            }
            crts[i] = CRT(p_p, crt_term);

        }
    } else if (r_cr==1){ //xi
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(p_l);
            for (int j = 0; j < p_l; j++){
                crt_term[j] = (2*r[i][j]+kd(j,i));
            }
            crts[i] = CRT(p_p, crt_term);
        }
    } else if (r_cr==2){ //ii
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(p_l);
            for (int j = 0; j < p_l; j++){
                crt_term[j] = (2*r[i][j]+(kd(j,i)*(power(2,p_rhoi+1))));
            }
            crts[i] = CRT(p_p, crt_term);
        }
    } else if (r_cr==3){ //o
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(p_l);
            for (int j = 0; j < p_l; j++){
                crt_term[j] = (2*r[i][j]+p_vert_s[i][j]);
            }
            crts[i] = CRT(p_p, crt_term);
        }
    } else { //op
        for(int i = 0; i < r_lenv; i++){
            std::vector<mpz_class> crt_term(p_l);
            for (int j = 0; j < p_l; j++){
                crt_term[j] = (2*r[i][j]+p_s[permute(j)][i]);
            }
            crts[i] = CRT(p_p, crt_term);
        }
    }

    
    for (int i = 0; i < r_lenv; i++){
        mpz_class chi_temp = floor_mod(r_chi.r_list[i],p_pi);

        deltas[i] = chi_temp+(E[i] * p_pi)-crts[i];
    }
    
    return deltas;
}

int Pk::permute(int j){
    if (j<6){
        return j+2;
    } else {
        return j-6;
    }
    
}

std::vector<mpz_class> Pk::make_front_U(PseudoRandomInts u_pri){
    std::vector<mpz_class> u(p_Theta);
    
    for (int i = 0; i < u_pri.r_len; i++){
        u[i] = u_pri.r_list[i];
    } //u draft
    
    int n = 0;
    for(int j = 0; j < p_l; j++){
        std::vector<int> s1indices;
        
        mpz_class xpj = floor_div(power(2, p_kap),p_p[j]); // i think its an int
        
        std::vector<mpz_class> su(p_Theta);
        for(int i = 0; i < p_Theta; i++){
            su[i] = (p_s[j][i] * u[i]);
            
            if (p_s[j][i] == 1){
                s1indices.push_back(i);
            }
        }
        
        mpz_class sumt = sum_array(su);
        sumt = floor_mod(sumt,power(2, p_kap+1));
        
        int v = n;
        n++;
        //change corresponding u
        su[v] = 0;
        mpz_class sumv = sum_array(su);
        
        mpz_class k1 = power(2, p_kap+1);
        mpz_class nu = k1 - sumv + xpj;
        
        while ((nu < 0) || (nu >= k1)){
            if (nu < 0){
                nu = nu+k1;
            } else {
                nu = nu-k1;
            }
        }
        
        u[v] = nu;

            
    }
    std::vector<mpz_class> front_u(p_l);
    
    for (int i = 0; i < p_l; i++){
        front_u[i] = u[i];
    }

    return front_u;
}
