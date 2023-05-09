all: mpi gcc

mpi:
	mpicc -o ./run_MPI.out ./projet_MPI.c -O3
	mpicc -o ./run_OMP.out ./projet_OMP.c -O3 -fopenmp

gcc:
	gcc -o ./run_SEQ.out ./projet_SEQ.c -O3
	gcc -o ./run_NAIVE.out ./projet_NAIVE.c -O3

run:
	$(info -----------Running with $(DIM) DIM------------)
	mpirun ./run_MPI.out $(DIM)
	mpirun -np 1 ./run_OMP.out $(DIM)
	./run_SEQ.out $(DIM)

run_with_naive:
	$(info -----------Running with $(DIM) DIM------------)
	mpirun ./run_MPI.out $(DIM)
	mpirun -np 1 ./run_OMP.out $(DIM)
	./run_SEQ.out $(DIM)
	./run_NAIVE.out $(DIM)
	
clean:
	rm -rf *.o
	rm -rf *.out