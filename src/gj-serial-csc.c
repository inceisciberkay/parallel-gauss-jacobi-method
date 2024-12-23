#include "utils.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 1.0e-10

void usage() {
  printf(
      "Usage: ./gaussJacobi-serial-csc <input_matrix_A_csc> <input_vector_B> "
      "<output_file>\n");
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Incorrect number of arguments.\n");
    usage();
    return 1;
  }

  const char *input_matrix_A_fname = argv[1];
  const char *input_vector_B_fname = argv[2];
  const char *output_vector_X_fname = argv[3];

  // read input matrix and vector
  int matrix_size, num_elements, *cols, *rows;
  double *elements; // matrix A
  double *b;        // vector b

  if (!read_input(input_matrix_A_fname, input_vector_B_fname, &matrix_size,
                  &num_elements, &cols, &rows, &elements, &b)) {
    fprintf(stderr, "Error while reading input.\n");
    return 1;
  }

  clock_t start = clock();

  // construct inital output vector X (calloc initializes values to 0)
  double *X_old = calloc(matrix_size, sizeof(double));

  // create a buffer to be used at middle steps (Y <- A*X_old && X_new <- Y + b)
  double *X_new =
      malloc(matrix_size * sizeof(double)); // Y will be refered as X_new

  // execute the algorithm
  double error = DBL_MAX;
  while (error > EPSILON) {
    // perform SpMV

    // in csc case, the values in X_new will not be computed sequentially.
    // Therefore a single accumulator (as in csr case) would not be sufficient
    // and the vector needs to be cleared at first to apply summations later.
    for (int i = 0; i < matrix_size; i++) {
      X_new[i] = 0;
    }

    for (int i = 0; i < matrix_size; i++) {
      for (int j = cols[i]; j < cols[i + 1]; j++) {
        X_new[rows[j]] += elements[j] * X_old[i];
      }
    }

    // remaining calculations
    double error_sum = 0;
    for (int i = 0; i < matrix_size; i++) {
      X_new[i] = X_new[i] + b[i];
      error_sum += pow(X_new[i] - X_old[i], 2);
    }

    error = sqrt(error_sum);

    // swap pointers
    double *temp = X_old;
    X_old = X_new;
    X_new = temp;
  }

  clock_t end = clock();
  double time = ((double)(end - start)) / CLOCKS_PER_SEC;
  printf("Time: %f seconds\n", time);

  // output result vector X to the file
  print_results(output_vector_X_fname, X_old, matrix_size);

  free(b);
  free(X_old);
  free(X_new);
  free(cols);
  free(rows);
  free(elements);

  return 0;
}
