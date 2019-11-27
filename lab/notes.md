### Following the Instructions

**Missing:**
+ Plots
+ Gauss-Seidel not Working in OMP

**Part 1**
1. OK
2. OK
3. No `eog` nor `gimp` installed. SCPing the image to local environment.
4. OK, all initial results saved in the `initial-results` directory.

**Part 2**
+ Error in `heat.h` line 24? Missing one algorithm in the comments?
+ I modify the `submit-omp.sh` script to take only powers of 2.
1. Jacobi
    + Parallelization of the Jacobi loop really does not do much?
        + Works a bit better for powers of 2
        + Too much overhead for really not much work maybe?
    + Why does it always do the same number of iterations?
    + Work a bit with different schedulers.
2. Red Black
    + It works way better
    + How many parallel regions do we have?
    + Could we merge both loops? Maybe nowait, one big loop?
        + Do as a second version (times)

**Part 3**
+ Initially confused by the fact that we lose track of the indices (given that we start allocating memory at 0, each process has its own vector).
1. Jacobi
    + This way, at any point, our row 0 and our row `myrows + 1` (in general), will belong to the neighbor. Hence, we need to receive from them their updated version, and send them ours.
    + It is important that we keep the same order to prevent deadlocks.
    + Further, we need to do some sanity checks in the master process to ensure he's not the only one.
    + Reduce + Bcast <=> Allreduce (google)
    + I know I could use Sendrcv but it was messy for me.
    + **Can't run with 8 processors! More processors than permitted**
    + I assumed Gauss-Seidel was not necessary (Tareador and stuff) although it should not be difficult.
2. Gauss:
    + Initially, well implemented but not speeding up.
    + **I feel that this is because we do an Allreduce with the global residual, PER ITERATION, hence we can not be at different iterations speeding up computation!!**
    + At the same time, if not all of them finish at the same time we can get stuck.
    + **Can't solve it :-(**

**Part 4**
+ I commented out the module extrae to load the software since it was giving problems, and we don't use it.
+ Changing the launcher aswell. Can't get it to work:
    - error while loading shared libraries: libcudart.so.4 
+ **Can't run!**
+ **Missing the shared memory thing!** -> Doubts
