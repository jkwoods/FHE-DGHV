//
//  utils.cpp
//  FHE
//
//  Created by Woods, Jess on 10/16/19.
//  work done at Oak Ridge National Lab
//

#include "utils.hpp"

#define w        (GMP_NUMB_BITS/2)
#define _BOTMASK    ((1ul << w)-1)
#define _TOPMASK    (~_BOTMASK)
#define BOT(x)        ((x) & _BOTMASK)
#define TOP(x)        ((x) >> w)
#define LIMB(z,i)    (((i)<((z)->_mp_size))?((z)->_mp_d[i]):(0L))
#define BOTL(z,i)    (BOT(LIMB(z,i)))
#define TOPL(z,i)    (TOP(LIMB(z,i)))
#define HLIMB(z,j)    ((j&1)?(TOPL(z,j>>1)):(BOTL(z,j>>1)))


mpq_class mod_2_f(mpq_class a){
    mpz_class num = a.get_num();
    mpz_class den = a.get_den()*2;
    mpz_class nearest_int = floor_div(num, den);

    mpq_class mod = a - nearest_int*2;

    
    if (mod < 0){
        while (mod <= 0){
            mod = mod + 2;
        }
    } else {
        while (mod >= 2){
            mod = mod - 2;
        }
    }
    return mod;
}

void print_vec(std::vector<int> p){
    for(int i = 0; i < p.size(); i++){
        std::cout << "i= " << i << " : " << p[i] << "\n";
    }
    std::cout << "\n";
}

void print_vec(std::vector<mpz_class> p){
    for(int i = 0; i < p.size(); i++){
        std::cout << "i= " << i << " : " << p[i] << "\n";
    }
    std::cout << "\n";
}

std::vector<int> to_binary(long a, int bits){ // [0,1,...,n] = [LSB, ... MSB]
    std::vector<int> result(bits);
    for(int i = 0; i < bits; i++){
        if (a == 0){
            break;
        }
        result[i] = a % 2;
        a = floor(a / 2);
    }
    return result;
}

mpz_class floor_mod(mpz_class a, int b){
    mpz_t r;
    mpz_init(r);

    mpz_t mb;
    mpz_init_set_ui(mb, 2);
    
    mpz_fdiv_r(r, a.get_mpz_t(), mb);
    
    mpz_class result = mpz_class(r);
    
    mpz_clear(r);
    return result;
}

mpz_class floor_mod(mpz_class a, mpz_class b){
    mpz_t r;
    mpz_init(r);
    mpz_fdiv_r(r, a.get_mpz_t(), b.get_mpz_t());
    
    mpz_class result = mpz_class(r);
    
    mpz_clear(r);
    return result;
}

mpz_class floor_div(mpz_class a, mpz_class b){
    mpz_t q;
    mpz_init(q);
    mpz_fdiv_q(q, a.get_mpz_t(), b.get_mpz_t());
    
    mpz_class result = mpz_class(q);
    
    mpz_clear(q);
    return result;
}

mpz_class modNear(mpz_class a, mpz_class b){ //convert to long/int?
    mpz_class q1 = (2*a+b);
    mpz_class q2 = (2*b);
    mpz_class quotientNear = floor_div(q1,q2); //q1 / q2;
    mpz_class mn = a - b*quotientNear;
    return mn;
}

mpz_class mul_inv(mpz_class a, mpz_class b){ //TODO - finish
    mpz_class b0 = b;
    mpz_class x0 = 0;
    mpz_class x1 = 1;
    if (b == 1){
        return 1;
    }
    while (a>1){
        mpz_class q = floor_div(a,b); //floor
        mpz_class temp = b;
        b = floor_mod(a,b);
        a = temp;
        
        mpz_class temp2 = x0;
        x0 = x1 - (q * x0);
        x1 = temp2;
    }
    if (x1 < 0){
        x1 += b0;
    }
    return x1;
}

mpz_class CRT(std::vector<mpz_class> n, std::vector<mpz_class> a){ //chinese remainder thm
    mpz_class prod = 1;
    for (int i = 0; i < n.size(); i++){
        prod = prod*n[i];
    }
    
    mpz_class sum = 0;
    for (int i = 0; i < n.size(); i++){
        mpz_class p = floor_div(prod,n[i]); //floor
        sum += a[i] * mul_inv(p, n[i]) * p;
    }

    mpz_class result = floor_mod(sum,prod);
    
    return result;
}

int kd(int i, int j){
  if (i == j){
    return 1;
  } else {
    return 0;
  }
}

mpz_class power(int base, int exp){
    mpz_t p;
    mpz_init(p);
    if (exp < 0){
         std::cout << "something is going wrong - negative exp in power function";
    }
    if (base < 0){ //handling for negative base (-2 mostly)
        base = base * (-1);
        mpz_ui_pow_ui(p, base, exp);

        mpz_class p_class(p);
        mpz_clear(p);
        
        return p_class*(-1);
    } else {
        mpz_ui_pow_ui(p, base, exp);

        mpz_class p_class(p);
        mpz_clear(p);
        
        return p_class;
    }

}

mpz_class random_prime_w(int ub, gmp_randstate_t rand_state){ //wierd range
    //generate mpz_t rand
    mpz_t p;
    mpz_init(p);
    mpz_urandomb(p, rand_state, ub-1); //0 - 2^(n-1)
    
    mpz_t ub_pow;
    mpz_init(ub_pow);
    mpz_ui_pow_ui(ub_pow, 2, ub-1);
    
    mpz_add(p, p, ub_pow);// + 2^(n-1)
    
    //check if prime
    while(mpz_probab_prime_p(p, 30) == 0){ //not prime
        mpz_urandomb(p, rand_state, ub-1); //0 - 2^(n-1)
        mpz_add(p, p, ub_pow);// + 2^(n-1)
    }
    
    //cast as mpz_class
    mpz_class p_class(p);
    
    mpz_clear(p);
    mpz_clear(ub_pow);
    
    return p_class;
}

mpz_class random_prime_f0(int ub, gmp_randstate_t rand_state){ //2^entry (only used one place)
    //generate mpz_t rand
    mpz_t p;
    mpz_init(p);
    mpz_urandomb(p, rand_state, ub);
    
    //check if prime
    while(mpz_probab_prime_p(p, 30) == 0){ //not prime
        mpz_urandomb(p, rand_state, ub);
    }
    
    //cast as mpz_class
    mpz_class p_class(p);
    mpz_clear(p);
    return p_class;
}




int random_choice(std::vector<int> sample){ //TODO test
    srand(time(0));
    int r = rand() % sample.size();
    return sample[r];
}

std::vector<int> random_sample(int range, int l){
    std::vector<int> sample(range);
    for(int i = 0; i < range; i++){
        sample[i] = i;
    }
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(sample.begin(), sample.end(), g);
    std::vector<int> cut_sample(l);
    for(int i = 0; i < l; i++){
        cut_sample[i] = sample[i];
    }
    return cut_sample;
}

mpz_class sum_array(std::vector<mpz_class> a){
    mpz_class suma = 0;
    for(int i = 0; i < a.size(); i++){
        suma = suma+a[i];
    }
    return suma;
}


std::vector<mpz_class> sum_binary(std::vector<mpz_class> a, std::vector<mpz_class> b){
    if (a.size() != b.size()){
        std::cout << "size issue in binary sum\n";
    }
    
    std::vector<mpz_class> c(a.size());
    c[0] = a[0]+b[0];
    //std::cout << "c[0]= " << c[0] << "\n";
    
    mpz_class carry = a[0]*b[0];
    //std::cout << "c[1] carry= " << carry << "\n";
    
    for (int i = 1; i < a.size()-1; i++){
        mpz_class carry2 = (a[i]+b[i])*carry+a[i]*b[i];
        mpz_class ci = a[i]+b[i]+carry;
        
        c[i] = (ci);
        //std::cout << "c[" << i << "]= " << c[i] << "\n";
        
        carry = carry2;
        //std::cout << "carry c[" << i+1 << "]= " << carry << "\n";
    }
    
    mpz_class ci = a.back()+b.back()+carry;
    c[a.size()-1] = ci;
    
    
    return c;
}

unsigned getGMP_NUMB_BITS()
{
  return GMP_NUMB_BITS;
}

unsigned long directScal(unsigned long kap, mpz_t cz, mpz_t yz)
{
    
  unsigned long nW=(kap+1)/(2*w), val=0, i;

  if(nW*w*2 != kap+1)
      return 0;

  for(i = 0; i < nW-1; i++) {
    val += BOTL(cz,i) * LIMB(yz,nW-1-i);
    val += (BOTL(cz,i) * TOPL(yz,nW-2-i)) >> w;
    val += TOPL(cz,i) * ((BOTL(yz,nW-1-i) << w) + TOPL(yz,nW-2-i));
    val += (TOPL(cz,i) * BOTL(yz,nW-2-i)) >> w;
  }

  val += BOTL(cz,nW-1) * LIMB(yz,0);
  val += TOPL(cz,nW-1) * (BOTL(yz,0) << w);

  return val;
}

