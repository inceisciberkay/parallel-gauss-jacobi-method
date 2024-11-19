import sys
from scipy import sparse
from scipy.sparse import random
import numpy as np

"""
Optionally you can alter the density by changing the default value of density.
"""


def gen_matrix_data(n, density=0.1):
    A = random(n, n, density=0.1, format="csr")
    A = A + (np.eye(n) + np.diagflat(A.sum(axis=1)))
    X = np.random.rand(n)
    B = np.matmul(A, X.T)
    D = np.diag(np.diag(A))
    R = A - D
    D_inv = np.linalg.inv(D)

    A_hat = np.matmul(-D_inv, R)
    b = np.matmul(D_inv, B.T)
    # X_test = np.linalg.solve(A, B.T)  # For checking the correction of X
    # print(X)
    # print(X_test.T)
    b = b.T.tolist()[0]  # Honestly I don't know why this is needed.
    b = np.array(b)  # If you see this don't judge me, if you ask why I can explain :D
    return A_hat, b, X


def print_matrix_to_file(matrix, name, dim="1D"):
    original_stdout = sys.stdout
    if dim == "1D":
        with open(name + ".txt", "w") as f:
            sys.stdout = f
            print(len(matrix))
            for i in matrix:
                print(i)

            sys.stdout = original_stdout

    elif dim == "2D_sparse":
        nnz = matrix.getnnz()
        indptr = matrix.indptr
        ind_len = len(indptr)
        indices = matrix.indices
        data = matrix.data
        with open(name + ".txt", "w") as f:
            sys.stdout = f
            print(ind_len - 1)
            print(nnz)
            for i in indptr:
                print(i)

            for i in indices:
                print(i)

            for i in data:
                print(i)

            sys.stdout = original_stdout
        pass


def main():
    if len(sys.argv) != 3:
        print("Proper usage: python DataGen.py [sampleName] [n]")
        print("[sampleName] is a custom file name for output files")
        print("n is the size of matrix A")
    else:
        A_hat, b, X = gen_matrix_data(int(sys.argv[2]))
        A_csr = sparse.csr_matrix(A_hat)
        A_csc = sparse.csc_matrix(A_hat)
        print_matrix_to_file(b, sys.argv[1] + "_matrix_b", "1D")
        print_matrix_to_file(X, sys.argv[1] + "_matrix_X", "1D")
        print_matrix_to_file(A_csr, sys.argv[1] + "_matrix_A_csr", "2D_sparse")
        print_matrix_to_file(A_csc, sys.argv[1] + "_matrix_A_csc", "2D_sparse")


if __name__ == "__main__":
    main()

