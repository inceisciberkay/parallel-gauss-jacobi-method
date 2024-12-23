all: gj-serial-csr gj-serial-csc gj-mpi-p1 gj-mpi-p2

gj-serial-csr: src/gj-serial-csr.c src/utils.c
	gcc src/gj-serial-csr.c src/utils.c -o gj-serial-csr -lm
	
gj-serial-csc: src/gj-serial-csc.c src/utils.c
	gcc src/gj-serial-csc.c src/utils.c -o gj-serial-csc -lm

gj-mpi-p1: src/gj-mpi-p1.c src/utils.c
	mpicc src/gj-mpi-p1.c src/utils.c -o gj-mpi-p1 -lm

gj-mpi-p2: src/gj-mpi-p2.c src/utils.c
	mpicc src/gj-mpi-p2.c src/utils.c -o gj-mpi-p2 -lm

clean:
	rm gj-serial-csr gj-serial-csc gj-mpi-p1 gj-mpi-p2
