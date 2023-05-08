all:
	mpicc ./projet.c

gcc:
	gcc ./projet_seq.c
	
clean:
	rm -rf *.o
	rm -rf *.out