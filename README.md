# "Fast" Matrix Multiplication

An attempt to implement various fast matrix multiplication algorithms in Julia, largely following [Markus Bläser's Theory of Computing Graduate Survey](https://theoryofcomputing.org/articles/gs005/gs005.pdf).

The word "fast" is in scare quotes because, while these algorithms _are_ asymptotically fast, in practice they are extremely slow, and these are far from optimized implementations.

I try to keep as much boilerplate as possible in `common.jl`.  The notebooks should contain only what's necessary to understand what makes one algorithm different from the next.