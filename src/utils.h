#ifndef UTIL_H_
#define UTIL_H_

// returns 1 if successful; 0 otherwise
int read_input(const char *A_fname, const char *b_fname, int *num_vertices,
               int *num_edges, int **vertices, int **edges, double **weights,
               double **b);

void print_results(const char *out_fname, double *output_vector_X,
                   int num_vertices);

// returns -1 if num is could not be found in arr
int find(int arr[], int size, int num);

#endif /* UTIL_H_ */
