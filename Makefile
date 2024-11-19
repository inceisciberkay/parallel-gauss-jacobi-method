all: gj-serial-csr gj-serial-csc gj-mpi-p1 gj-mpi-p2

gj-serial-csr: src/gj-serial-csr.c src/util.c
	gcc src/gj-serial-csr.c src/util.c -o gj-serial-csr -lm
	
gj-serial-csc: src/gj-serial-csc.c src/util.c
	gcc src/gj-serial-csc.c src/util.c -o gj-serial-csc -lm

gj-mpi-p1: src/gj-mpi-p1.c src/util.c
	mpicc src/gj-mpi-p1.c src/util.c -o gj-mpi-p1 -lm

gj-mpi-p2: src/gj-mpi-p2.c src/util.c
	mpicc src/gj-mpi-p2.c src/util.c -o gj-mpi-p2 -lm

clean:
	rm gj-serial-csr gj-serial-csc gj-mpi-p1 gj-mpi-p2 out.txt
