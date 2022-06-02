# Naive-Gaussian-Elimination
A program that computes solutions for linear systems. 
Takes as input a file which contains data for a linear system in the following format:

n
a11 a12 ... a1n
a21 a22 ... a2n
...
an1 an2 ... ann
b1   b2   ... bn

The user is able to modify the programs behavior with optional flag --spp, in which case the program will use Scaled Partial Pivoting to produce the solution. For example, for a system placed in file sys1.lin, the user could run:

> gaussian sys1.lin

or,

> gaussian --spp sys1.lin

In the first case the program will use Naive-Gaussian Elimination, and in the second it will use SPP. In both cases, the solution will be placed in file sys1.sol.
