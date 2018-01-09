# Cholesky decomposition (Romanian language)



    Arhitecturi si Prelucrari Paralele (C1)

    Comenzi de compilare si de rulare pentru fiecare versiune (Mac OS)

    Badila Gabriel Alin
    Craciun Alexandru Sever
    342 C1


    =========================================================================================

    Varianta Seriala
        compilare:  gcc-7 -o serial serialCholesky.c
        rulare:     ./serial <in_file>

    Varianta OpenMP
        compilare:  g++-7 -o openmp openmpCholesky.cpp -fopenmp
        rulare:     ./openmp <in_file>

    Varianta MPI
        compilare:  mpicc -o mpi mpiCholesky.c
        rulare:     mpirun --oversubscribe -np <num_procs> ./mpi <in_file>

    Varianta Pthreads
        compilare:  gcc-7 -o pthreads pthreadsCholesky.c Utils/barrier.c -lpthread -lm
        rulare:     ./pthreads <in_file>

    =========================================================================================
