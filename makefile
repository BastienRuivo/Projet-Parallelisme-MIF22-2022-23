all: mpicc gcc

mpi:
	mpicc -o  ./run ./projet.c -O3

gcc:
	gcc -o  ./run_seq ./projet_seq.c -O3
	
clean:
	rm -rf *.o
	rm -rf *.out