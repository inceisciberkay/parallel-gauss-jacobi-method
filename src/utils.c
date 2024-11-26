#include "utils.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int read_input_matrix_A_file(const char *file_name, int *num_vertices,
                             int *num_edges, int **vertices, int **edges,
                             double **weights) {
  printf("Reading matrix A file...\n");
  FILE *input_file = fopen(file_name, "r");
  if (input_file) {
    int success;

    success = fscanf(input_file, "%d", &*num_vertices);
    if (!success) {
      printf("Bad File format!\n");
      return 0;
    }
    success = fscanf(input_file, "%d", &*num_edges);
    if (!success) {
      printf("Bad File format!\n");
      return 0;
    }

    // topologicalSize includes the edge weights
    int topologicalSize = (*num_vertices + 1 + (*num_edges * 2)) * sizeof(int);

    printf("num_vertices = %d, num_edges = %d\n", *num_vertices, *num_edges);

    printf("Graph data footprint=%.3f KB ~ %.3f MB ~ %.3f GB\n",
           topologicalSize / 1024.0, topologicalSize / (1024.0 * 1024),
           topologicalSize / (1024.0 * 1024 * 1024));

    *vertices = (int *)malloc((*num_vertices + 1) * sizeof(int));
    *edges = (int *)malloc(*num_edges * sizeof(int));
    *weights = (double *)malloc(*num_edges * sizeof(double));

    // read vertices
    for (int i = 0; i < *num_vertices + 1; i++) {
      success = fscanf(input_file, "%d", &((*vertices)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    // read edges
    for (int i = 0; i < *num_edges; i++) {
      success = fscanf(input_file, "%d", &((*edges)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    // read weights
    for (int i = 0; i < *num_edges; i++) {
      success = fscanf(input_file, "%lf", &((*weights)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    fclose(input_file);
    return 1;
  }
  printf("Could not open the file!\n");
  return 0;
}

int read_input_vector_b_file(const char *fname, double **b) {
  FILE *input_file_b = fopen(fname, "r");

  if (!input_file_b) {
    fprintf(stderr,
            "File for input vector b could not be opened for reading.\n");
    return 0;
  }

  int num_vertices;
  if (!fscanf(input_file_b, "%d", &num_vertices)) {
    fprintf(stderr, "First line of input file for vector b should be the "
                    "number of vertices.\n");
    return 0;
  };

  *b = malloc(num_vertices * sizeof(double));

  for (int i = 0; i < num_vertices; i++) {
    if (!fscanf(input_file_b, "%lf", *b + i)) {
      fprintf(stderr, "Bad file format.\n");
      return 0;
    }
  }

  fclose(input_file_b);
  return 1;
}

int read_input(const char *A_fname, const char *b_fname, int *num_vertices,
               int *num_edges, int **vertices, int **edges, double **weights,
               double **b) {

  // read input matrix A
  if (!read_input_matrix_A_file(A_fname, num_vertices, num_edges, vertices,
                                edges, weights)) {
    fprintf(stderr, "Input matrix A file could not be read.\n");
    return 0;
  }

  // read input vector b
  if (!read_input_vector_b_file(b_fname, b)) {
    fprintf(stderr, "Input vector b file could not be read.\n");
    return 0;
  }

  return 1;
}

void print_results(const char *out_fname, double *output_vector_X,
                   int num_vertices) {
  FILE *output_file = fopen(out_fname, "w");
  if (output_file) {
    // output the number of vertices
    fprintf(output_file, "%d\n", num_vertices);
    // Output the final values
    for (int i = 0; i < num_vertices; i++) {
      fprintf(output_file, "%.17lf\n", output_vector_X[i]);
    }
    fclose(output_file);
  } else {
    printf("Could not open the file!\n");
  }
}

int find(int arr[], int size, int num) {
  for (int i = 0; i < size; i++) {
    if (arr[i] == num)
      return i;
  }
  return -1;
}
