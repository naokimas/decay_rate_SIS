# decay_rate_SIS
Decay rate and its bound for the susceptible-infectious-susceptible (SIS) model on networks

This code produces 
- the decay rate of the number of infectious nodes in the SIS model by direct stochastic numerical simulations of the SIS model on networks, using the Gillespie's direct algorithm, and
- the theoretical bounds of the decay rate.

If you use this code, please cite

[Naoki Masuda, Victor M. Preciado, Masaki Ogura.  
Analysis of the susceptible-infected-susceptible epidemic dynamics in networks via the non-backtracking matrix.  
IMA Journal of Applied Mathematics, 85, 214-230 (2020).](https://doi.org/10.1093/imamat/hxaa003)

This paper is open access.

Usage (i)

```
g++ sis-net-decay.cc  
a.out input-file-name
```

- You need a C++ compiler.
- We provide four edge list files, any of which can be used as input-file-name.
- The input file has to follow a certain format, which is described in the preample of sis-net-decay.cc
- One can alternatively set the data_type as an integer in place of input-file-name when running a.out. See sis-net-decay.cc for details.

Usage (ii)
