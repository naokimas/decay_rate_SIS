# Decay rate and its theoretical bounds for the susceptible-infectious-susceptible (SIS) model on networks

- `sis-net-decay.cc` is a C/C++ code and produces the decay rate of the number of infectious nodes in the SIS model by direct stochastic numerical simulations of the SIS model on networks, using the Gillespie's direct algorithm.
- `siso.m` is a MATLAB code and produces the theoretical bounds of the decay rate.

If you use this code, please cite

[Naoki Masuda, Victor M. Preciado, Masaki Ogura.  
Analysis of the susceptible-infected-susceptible epidemic dynamics in networks via the non-backtracking matrix.  
IMA Journal of Applied Mathematics, 85, 214-230 (2020).](https://doi.org/10.1093/imamat/hxaa003)

This paper is open access.

## Calculate the decay rate by direct stochastic numerical simulations of the SIS model on networks

```
g++ sis-net-decay.cc  
a.out input-file-name
```

- You need a C++ compiler.
- We provide four edge list files (with extension .mat, but they are not MATLAB files; they are text files) in this folder. Any of these files can be specified as input-file-name. These four networks are used in the paper.
- The input file has to follow a certain format, which is described in the comment lines in the beginning of sis-net-decay.cc
- One can alternatively set the data_type as an integer in place of input-file-name when running a.out. See the code for details.

## Calculate theoretical lower bounds of the decay rate of the SIS model on networks

Run `siso.m` on MATLAB without arguments.

- It reads an edge list file and outputs a table of theoretical lower bounds of the decay rate for infection rate values.

