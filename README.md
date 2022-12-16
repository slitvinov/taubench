<h2>TauBench</h2>

This is an unstructured grid benchmark. The benchmark mimics the
run-time performance of the Tau solver (a Navier Stokes solver which
has been developed at the DLR in Germany). Tau is a three-dimensional
parallel hybrid multigrid solver, which uses a finite volume scheme in
order to solve the Reynolds-averaged Navier-Stokes equations.

<h2>Install</h2>
<pre>
$ make
</pre>


To change the compiler

<pre>
$ make MPICC=mpicc.mpich
</pre>

<h2>Run</h2>
Tau supports both vector and cache colored grids, you might also need
to adapt the compiler directives in
<a href="./nodep.h">nodep.h</a>
and
<a href="./expand.h">expand.h</a>.
The benchmark itself can be used with two flags, the gridsize per
process and the number of pseudo steps:

<pre>
$ mpirun -np 2 ./taubench -n 100000 -s 10
This is TauBench.
Evaluating kernels - please be patient.
..........

        - kernel_1_0 :      0.934 secs -   1506.596 mflops
        - kernel_1_1 :      0.349 secs -    625.287 mflops
        - kernel_2_1 :      0.619 secs -   1603.554 mflops
        - kernel_2_2 :      0.809 secs -   1508.574 mflops
        - kernel_2_3 :      0.310 secs -    794.652 mflops
        - kernel_2_4 :      0.472 secs -   1174.215 mflops
        - kernel_3_0 :      3.048 secs -    898.502 mflops

               total :      6.703 secs -   1888.867 mflops

points     :     100000
steps      :         10
procs      :          2

comp       :      6.351 secs
comm       :      0.352 secs
comm ratio :      0.055
</pre>
