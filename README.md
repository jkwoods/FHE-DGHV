# FHE-DGHV
Parallel C++ Implementation of the Fully Homomorphic Batch DGHV scheme with compressed public key <br />
<br />

**Notes** <br />
This was made for the purpose of fun/learning/proof of concept, and should NOT be used for any secure application.<br />

Please feel free to use/modify code. Everything is open-source.

**Pre-Req**: <br />
GMP, The GNU Multiple Precision Arithmetic Library https://gmplib.org/
CUDA GPU offloading is available, but not required.
OpenMP parallelization is available, but not required.

**Using Library** <br />
Keys can be created with preset parameters: <br />
```Pk example_pk = Pk::make_key(int security_level);``` <br />
Security levels toy, small, medium and large, are inputs 0-3, respectively. <br />
You can also set your own parameters with the Pk constructor: <br />
```Pk example_pk = Pk(int lam, int rho, int eta, int gam, int Theta, int alpha, int tau, int l);```
```example_pk.assert_parameter_correctness();```

Data can be encoded with the "Encoding" class:<br />
```std::vector<int> data = {0,1,0,0,0,1,1,0,1 ...};```
```Encoding example_encoding = Encoding(Pk publicKey, std::vector<int> data);``` <br />
All encoded data should be treated as a vector of l (seperate) bits.

Addition and multiplication have been overloaded. They are done bitwise. <br />

Negation (of each individual "slot") is available:<br />
```example_encoding.neg();``` <br />
As is selection (as a static method):<br />
```Encoding example_selected = Encoding::selector(std::vector<int> s, Encoding a, Encoding b); //if s=1: a, else b```

Recoding happens automatically after enough operations, but can also be done explicitly:<br />
```example_encoding.recode();```


**References/Acknowledgements** <br />
Sequential version based on these papers: <br />
https://eprint.iacr.org/2013/036.pdf <br />
https://eprint.iacr.org/2011/441.pdf <br />
https://eprint.iacr.org/2011/440.pdf <br />
https://eprint.iacr.org/2009/616.pdf <br />

And loosely on this (sequential, non-batched) code:
https://github.com/coron/fhe <br />

GMP Library used to handle big integers. <br />

This work was done at Oak Ridge National Lab. This work was supported in part by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists (WDTS) under the Science Undergraduate Laboratory Internship program.	

Contact: jesskwoods (at) gmail.com


