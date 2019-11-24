# FHE-DGHV
Parallel C++ Implementation of the Fully Homomorphic Batch DGHV scheme with compressed public key <br />
<br />

**Notes** <br />
This was made for the purpose of fun/learning/proof of concept, and should NOT be used for any secure application.<br />

Please feel free to use/modify code. Everything is open-source.

**Using Library** <br />
Keys can be made by creating a "Pk" object:


Data can be encoded with the "Encoding" class:
```Encoding example = Encoding(Pk publicKey, std::vector<int> {1,0,0,1,0,...});```

Addition and multiplication have been overloaded.


**References/Acknowledgements** <br />
Sequential version based on these papers: <br />
https://eprint.iacr.org/2013/036.pdf <br />
https://eprint.iacr.org/2011/441.pdf <br />
https://eprint.iacr.org/2011/440.pdf <br />
https://eprint.iacr.org/2009/616.pdf <br />

And loosely on this (sequential, non-batched) code:
https://github.com/coron/fhe <br />

GMP Library (https://gmplib.org/) used to handle big integers. <br />

This work was done at Oak Ridge National Lab. This work was supported in part by the U.S. Department of Energy, Office of Science, Office of Workforce Development for Teachers and Scientists (WDTS) under the Science Undergraduate Laboratory Internship program.	

Contact: jesskwoods (at) gmail.com


