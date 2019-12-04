# FHE-DGHV
Parallel C++ Implementation of the Fully Homomorphic Batch DGHV scheme with compressed public key <br />
<br />

**Notes** <br />
This was made for the purpose of fun/learning/proof of concept, and should NOT be used for any secure application.<br />

Please feel free to use/modify code. Everything is open-source.

**Pre-Req**: <br />
GMP, The GNU Multiple Precision Arithmetic Library https://gmplib.org/ <br />

```
./configure --enable-cxx   //more complex options for specific processors/ISAs are available through the GMP manual
make
make check 
make install
```


CUDA GPU offloading is available, but not required. <br />
OpenMP parallelization is available, but not required. <br />

FHE-DGHV Install:
```
git clone https://github.com/jkwoods/FHE-DGHV
cd FHE-DGHV
make    //TODO: set gmp library path
make clean
```

**Using Library** <br />
Keys can be created with preset parameters: <br />
```
Pk example_pk = Pk::make_key(int security_level);
```
Security levels toy, small, medium and large, are inputs 0-3, respectively. <br />
You can also set your own parameters with the Pk constructor: <br />
```
Pk example_pk = Pk(int lam, int rho, int eta, int gam, int Theta, int alpha, int tau, int l);
example_pk.assert_parameter_correctness();
```

Data can be encoded with the "Encoding" class:<br />
```
std::vector<int> data = {0,1,0,0,0,1,1,0,1 ...};
Encoding example_encoding = Encoding(Pk publicKey, std::vector<int> data);
```
All encoded data should be treated as a vector of l (seperate) bits.

Addition and multiplication have been overloaded. They are done bitwise. <br />

Negation (of each individual "slot") is available:<br />
```
example_encoding.neg();
```
As is selection (as a static method):<br />
```
Encoding example_selected = Encoding::selector(std::vector<int> s, Encoding a, Encoding b); //if s=1: a, else b
```

Recoding happens automatically after enough operations, but can also be done explicitly:<br />
```
example_encoding.recode();
```


**References/Acknowledgements** <br />
Sequential version based on these papers (cited below): <br />
https://eprint.iacr.org/2013/036.pdf <br />
https://eprint.iacr.org/2011/441.pdf <br />
https://eprint.iacr.org/2011/440.pdf <br />
https://eprint.iacr.org/2009/616.pdf <br />

And loosely on this (sequential, non-batched) code:
https://github.com/coron/fhe <br />

GMP Library used to handle big integers. <br />

This work was done at Oak Ridge National Lab. This work was supported in part by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists (WDTS) under the Science Undergraduate Laboratory Internship program.	

Contact: jesskwoods (at) gmail.com

**Bibliography** <br />
Craig Gentry et al. Fully homomorphic encryption using ideal lattices. InStoc, volume 9,pages 169–178, 2009.

Craig  Gentry  and  Shai  Halevi.   Implementing  gentrys  fully-homomorphic  encryptionscheme.  InAnnual international conference on the theory and applications of crypto-graphic techniques, pages 129–148. Springer, 2011.

Marten  Van  Dijk,  Craig  Gentry,  Shai  Halevi,  and  Vinod  Vaikuntanathan.   Fully  ho-momorphic  encryption  over  the  integers.   InAnnual International Conference on theTheory and Applications of Cryptographic Techniques, pages 24–43. Springer, 2010.

Jean-S ́ebastien Coron, David Naccache, and Mehdi Tibouchi.  Public key compressionand modulus switching for fully homomorphic encryption over the integers.  InAnnualInternational Conference on the Theory and Applications of Cryptographic Techniques,pages 446–464. Springer, 2012.

Nigel P Smart and Frederik Vercauteren. Fully homomorphic simd operations.Designs,codes and cryptography, 71(1):57–81, 2014.

Jung  Hee  Cheon,  Jean-S ́ebastien  Coron,  Jinsu  Kim,  Moon  Sung  Lee,  Tancrede  Le-point,  Mehdi  Tibouchi,  and  Aaram  Yun.   Batch  fully  homomorphic  encryption  overthe  integers.   InAnnual International Conference on the Theory and Applications ofCryptographic Techniques, pages 315–335. Springer, 2013.

Jean-S ́ebastien Coron, Avradip Mandal, David Naccache, and Mehdi Tibouchi.  Fullyhomomorphic encryption over the integers with shorter public keys.  InAnnual Cryp-tology Conference, pages 487–504. Springer, 2011.

Python Core Team.Python: A dynamic, open source programming language.PythonSoftware Foundation, 2015.

Implementationofthedghvfullyhomomorphicencryptionscheme.https://github.com/coron/fhe, 2012.

Dask Development Team.Dask: Library for dynamic task scheduling, 2016.

Bjarne Stroustrup.The C++ programming language.  Pearson Education, 2013.

Torbjrn  Granlund  and  the  GMP  development  team.GNU MP: The GNU MultiplePrecision Arithmetic Library, 5.0.5 edition, 2012.11

OpenMP Architecture Review Board.  OpenMP application program interface version3.0, May 2008.

Simon Fau, Renaud Sirdey, Caroline Fontaine, Carlos Aguilar-Melchor, and Guy Gog-niat. Towards practical program execution over fully homomorphic encryption schemes.In2013 Eighth International Conference on P2P, Parallel, Grid, Cloud and InternetComputing, pages 284–290. IEEE, 2013.

Craig  Gentry,  Shai  Halevi,  and  Nigel  P  Smart.   Homomorphic  evaluation  of  the  aescircuit.  InAnnual Cryptology Conference, pages 850–867. Springer, 2012.

Jungwon Kim and Jeffrey S Vetter.  Implementing efficient data compression and en-cryption  in  a  persistent  key-value  store  for  hpc.The International Journal of HighPerformance Computing Applications, 33(6):1098–1112, 2019.

Wei Wang, Yin Hu, Lianmu Chen, Xinming Huang, and Berk Sunar. Accelerating fullyhomomorphic  encryption  using  gpu.   In2012 IEEE conference on high performanceextreme computing, pages 1–5. IEEE, 2012.

Wei Wang, Yin Hu, Lianmu Chen, Xinming Huang, and Berk Sunar. Exploring the fea-sibility of fully homomorphic encryption.IEEE Transactions on Computers, 64(3):698–706, 2013.

J. Nickolls.  Scalable parallel programming with cuda introduction.  In 2008 IEEE HotChips 20 Symposium (HCS), pages 1–9, Aug 2008.

Jeff  Bezanson,  Alan  Edelman,  Stefan  Karpinski,  and  Viral  B  Shah.   Julia:   A  freshapproach to numerical computing.SIAM review, 59(1):65–98, 2017.


