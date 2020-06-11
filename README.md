## TauBench

This is the alpha release of a unstructured grid benchmark. TauBench
is a pseudo benchmark. The respective kernels are derived from Tau --
a Navier Stokes solver which has been developed at the German
aerospace research DLR in germany. The original flow solver is a
three-dimensional parallel hybrid multigrid solver, which uses a
finite volume scheme in order to solve the Reynolds-averaged
Navier-Stokes equations.  The benchmark itself doesn't do anything
useful in terms of producing meaningful results - rather it is meant
to mimic the run-time performance of the Tau solver.

## INSTALL
~~~
make
~~~

Tau supports both vector and cache colored grids, you might also need
to adapt the compiler directives in nodep.h and expand.h.

The benchmark itself can be used with two flags, the gridsize per
process and the number of pseudo steps:

~~~
mpirun -np 2 ./taubench -n 100000 -s 10
~~~
