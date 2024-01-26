<h2>TauBench</h2>
This is an unstructured grid benchmark. The benchmark mimics the
run-time performance of the Tau solver (a Navier Stokes solver which
has been developed at the DLR in Germany). Tau is a three-dimensional
parallel hybrid multigrid solver, which solves the Reynolds-averaged
Navier-Stokes equations by finite volume scheme.

<h2>Install</h2>
<pre>
$ make
mpicc -g -O2 flux.c -c
mpicc -g -O2 main.c -c
mpicc -g -O2 lim.c -c
mpicc -g -O2 smooth.c -c
mpicc flux.o main.o lim.o smooth.o   -lm -o taubench
</pre>

To change the compiler and compiler flags

<pre>
$ make 'MPICC = mpicc.mpich' 'CFLAGS = -Ofast -Wall'
mpicc.mpich -Ofast -Wall flux.c -c
mpicc.mpich -Ofast -Wall main.c -c
mpicc.mpich -Ofast -Wall lim.c -c
mpicc.mpich -Ofast -Wall smooth.c -c
mpicc.mpich flux.o main.o lim.o smooth.o   -lm -o taubench
</pre>

<h2>Run</h2>
Tau supports both vector and cache colored grids, you might also need
to adapt the compiler directives in
<a href="./nodep.inc">nodep.inc</a>
and
<a href="./expand.inc">expand.inc</a>.
The benchmark itself can be used with two flags, the gridsize per
process and the number of pseudo steps:

<pre>
$ mpiexec -n 2 ./taubench -n 100000 -s 10
This is TauBench.
Evaluating kernels - please be patient.
..........

        - kernel_1_0 :      0.315 secs -   4470.269 mflops
        - kernel_1_1 :      0.131 secs -   1666.135 mflops
        - kernel_2_1 :      0.210 secs -   4733.000 mflops
        - kernel_2_2 :      0.186 secs -   6560.788 mflops
        - kernel_2_3 :      0.089 secs -   2759.613 mflops
        - kernel_2_4 :      0.069 secs -   7992.762 mflops
        - kernel_3_0 :      0.344 secs -   7952.712 mflops

               total :      1.311 secs -   9658.964 mflops

points     :     100000
steps      :         10
procs      :          2

comp       :      1.302 secs
comm       :      0.009 secs
comm ratio :      0.007
</pre>
