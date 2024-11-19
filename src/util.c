#include "util.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int read_input_matrix_A_file(char *file_name, int *numVertices, int *numEdges,
                             int **vertices, int **edges, double **weights) {
  printf("Reading matrix A file...\n");
  FILE *inputFile = fopen(file_name, "r");
  if (inputFile) {
    int success;

    success = fscanf(inputFile, "%d", &*numVertices);
    if (!success) {
      printf("Bad File format!\n");
      return 0;
    }
    success = fscanf(inputFile, "%d", &*numEdges);
    if (!success) {
      printf("Bad File format!\n");
      return 0;
    }

    // topologicalSize includes the edge weights
    int topologicalSize = (*numVertices + 1 + (*numEdges * 2)) * sizeof(int);

    printf("numVertices = %d, numEdges = %d\n", *numVertices, *numEdges);

    printf("Graph data footprint=%.3f KB ~ %.3f MB ~ %.3f GB\n",
           topologicalSize / 1024.0, topologicalSize / (1024.0 * 1024),
           topologicalSize / (1024.0 * 1024 * 1024));

    *vertices = (int *)malloc((*numVertices + 1) * sizeof(int));
    *edges = (int *)malloc(*numEdges * sizeof(int));
    *weights = (double *)malloc(*numEdges * sizeof(double));

    // read vertices
    for (int i = 0; i < *numVertices + 1; i++) {
      success = fscanf(inputFile, "%d", &((*vertices)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    // read edges
    for (int i = 0; i < *numEdges; i++) {
      success = fscanf(inputFile, "%d", &((*edges)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    // read weights
    for (int i = 0; i < *numEdges; i++) {
      success = fscanf(inputFile, "%lf", &((*weights)[i]));
      if (success == EOF || success == 0) {
        printf("Bad File format!\n");
        return 0;
      }
    }

    fclose(inputFile);
    return 1;
  }
  printf("Could not open the file!\n");
  return 0;
}

int read_input_vector_b_file(char *filename, double **b) {
  FILE *input_file_b = fopen(filename, "r");

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

void printResults(char *fileName, double *output_vector_X, int numVertices) {
  FILE *output_file = fopen(fileName, "w");
  if (output_file) {
    // output the number of vertices
    fprintf(output_file, "%d\n", numVertices);
    // Output the final values
    for (int i = 0; i < numVertices; i++) {
      fprintf(output_file, "%.17lf\n", output_vector_X[i]);
    }
    fclose(output_file);
  } else {
    printf("Could not open the file!\n");
  }
}
