# Mutually Unbiased Bases

Construction of mutually unbiased bases (MUBs) in various languages.
Each language works independently of the others.
Any contribution is welcome!

## Languages supported

### Julia
The code has been improved and is now part of the package [Ket](https://github.com/araujoms/Ket.jl), which should be preferred in all cases, both for the increased quality of the implementation of these MUBs functions, but also for the ease of installation and the numerous other functionalities, for instance, the construction of symmetric, informationnally complete, positive operator-valued measures (SIC-POVMs).

Dependencies: [Combinatorics](https://juliamath.github.io/Combinatorics.jl/dev/), [Nemo](https://nemocas.github.io/Nemo.jl/stable/) (there was a breaking update at some point, manually select the correct line if you run into problems), [Primes](https://juliamath.github.io/Primes.jl/stable/)

### Mathematica
Dependency: [FiniteFields](https://reference.wolfram.com/language/FiniteFields/guide/FiniteFieldsPackage.html)

### Matlab
Dependency: [gf](https://de.mathworks.com/matlabcentral/fileexchange/32872-a-toolbox-for-simple-finite-field-operation)

## References
Stephen Brierley, Stefan Weigert, Ingemar Bengtsson, *All mutually unbiased bases in dimensions two to five*, [arXiv:0907.4097](https://arxiv.org/abs/0907.4097), [Quantum Information & Computation **10**, 0803–0820 (2010)](https://doi.org/10.26421/QIC10.9-10-6).\
Thomas Durt, Berthold-Georg Englert, Ingemar Bengtsson, and Karol Życzkowski, *On mutually unbiased bases*, [arXiv:1004.3348](https://arxiv.org/abs/1004.3348), [International Journal of Quantum Information **8**, 535–640 (2010)](https://doi.org/10.1142/S0219749910006502).\
William K. Kantor, *MUBs inequivalence and affine planes*, [arXiv:1104.3370](https://arxiv.org/abs/1104.3370), [Journal of Mathematical Physics **53**, 032204 (2012)](https://doi.org/10.1063/1.3690050).
