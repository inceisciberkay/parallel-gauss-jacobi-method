#include "util.h"
#include <float.h>
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define EPSILON 1.0e-10

void usage() {
  printf("Usage: ./gaussJacobi-mpi-p1 <input_matrix_A_csr> <input_vector_B> "
         "<output_file>\n");
}

// returns -1 if not found
int find(int arr[], int size, int num) {
  for (int i = 0; i < size; i++) {
    if (arr[i] == num)
      return i;
  }
  return -1;
}

int main(int argc, char *argv[]) {
  if (argc != 4) {
    fprintf(stderr, "Incorrect number of arguments.\n");
    usage();
    return 1;
  }

  int myid, numprocs;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  int matrix_size, num_values, *rows, *columns;
  double *values;
  double *b;

  if (myid == 0) { // master process
    // reading input files
    char *input_matrix_A_filename = argv[1];
    char *input_vector_B_filename = argv[2];

    // read input matrix A
    if (!read_input_matrix_A_file(input_matrix_A_filename, &matrix_size,
                                  &num_values, &rows, &columns, &values)) {
      fprintf(stderr, "Input matrix A file could not be read.\n");
      return 1;
    }

    // read input vector b
    if (!read_input_vector_b_file(input_vector_B_filename, &b)) {
      fprintf(stderr, "Input vector b file could not be read.\n");
      return 1;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
  double start = MPI_Wtime();

  // master process broadcasts the matrix size so that each processor can
  // compute its recvsize
  MPI_Bcast(&matrix_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

  const int LOCAL_ROW_COUNT = matrix_size / numprocs;

  /* distribute input matrix A and input vector b to child processes */

  // partitioning rows, columns and values (uniform block row-wise
  // partitioning)
  int *send_counts_rows, *displs_rows, *send_counts_columns, *displs_columns,
      *send_counts_b, *displs_b;
  if (myid == 0) {
    send_counts_rows = malloc(numprocs * sizeof(int));
    displs_rows = malloc(numprocs * sizeof(int));

    send_counts_columns = malloc(numprocs * sizeof(int));
    displs_columns = malloc(numprocs * sizeof(int));

    send_counts_b = malloc(numprocs * sizeof(int));
    displs_b = malloc(numprocs * sizeof(int));

    for (int i = 0; i < numprocs; i++) {
      send_counts_rows[i] = LOCAL_ROW_COUNT + 1;
      displs_rows[i] = i * LOCAL_ROW_COUNT;

      send_counts_columns[i] =
          rows[(i + 1) * (LOCAL_ROW_COUNT)] - rows[i * (LOCAL_ROW_COUNT)];
      displs_columns[i] =
          i == 0 ? 0 : displs_columns[i - 1] + send_counts_columns[i - 1];

      send_counts_b[i] = LOCAL_ROW_COUNT;
      displs_b[i] = i * LOCAL_ROW_COUNT;
    }
  }

  int rcvsize_rows = LOCAL_ROW_COUNT + 1;
  int *local_rows = malloc(rcvsize_rows * sizeof(int));
  MPI_Scatterv(rows, send_counts_rows, displs_rows, MPI_INT, local_rows,
               rcvsize_rows, MPI_INT, 0, MPI_COMM_WORLD);

  int rcvsize_columns = local_rows[rcvsize_rows - 1] - local_rows[0];
  int *local_columns = malloc(rcvsize_columns * sizeof(int));
  double *local_values = malloc(rcvsize_columns * sizeof(double));
  MPI_Scatterv(columns, send_counts_columns, displs_columns, MPI_INT,
               local_columns, rcvsize_columns, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(values, send_counts_columns, displs_columns, MPI_DOUBLE,
               local_values, rcvsize_columns, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // partitioning input vector b
  int rcvsize_b = LOCAL_ROW_COUNT;
  double *local_b = malloc(rcvsize_b * sizeof(double));
  MPI_Scatterv(b, send_counts_b, displs_b, MPI_DOUBLE, local_b, rcvsize_b,
               MPI_DOUBLE, 0, MPI_COMM_WORLD);

  if (myid == 0) {
    free(rows);
    free(columns);
    free(values);
    free(b);
    free(send_counts_rows);
    free(displs_rows);
    free(send_counts_columns);
    free(displs_columns);
    free(send_counts_b);
    free(displs_b);
  }

  // constructing receive list
  int myrecv_list_size = 0;
  int *myrecv_list = NULL;
  int *myrecv_counts = NULL;
  int **myrecv_indexes = NULL;

  for (int i = 0; i < LOCAL_ROW_COUNT; i++) {
    int displ = local_rows[0];
    for (int k = local_rows[i]; k < local_rows[i + 1]; k++) {
      int j = local_columns[k - displ];
      int jpid = j / (LOCAL_ROW_COUNT); // mapping
      if (jpid != myid) {
        int index_of_jpid_in_myrecv_list =
            find(myrecv_list, myrecv_list_size, jpid);

        if (index_of_jpid_in_myrecv_list == -1) {
          myrecv_list_size++;

          myrecv_list = realloc(myrecv_list, myrecv_list_size * sizeof(int));
          myrecv_counts =
              realloc(myrecv_counts, myrecv_list_size * sizeof(int));
          myrecv_indexes =
              realloc(myrecv_indexes, myrecv_list_size * sizeof(int *));

          index_of_jpid_in_myrecv_list = myrecv_list_size - 1;
          myrecv_list[index_of_jpid_in_myrecv_list] = jpid;
          myrecv_counts[index_of_jpid_in_myrecv_list] = 0;
          myrecv_indexes[index_of_jpid_in_myrecv_list] = NULL;
        }

        // check if j is already included in myrecv_indexes
        if (find(myrecv_indexes[index_of_jpid_in_myrecv_list],
                 myrecv_counts[index_of_jpid_in_myrecv_list], j) == -1) {
          // if not included
          myrecv_counts[index_of_jpid_in_myrecv_list]++;
          myrecv_indexes[index_of_jpid_in_myrecv_list] = realloc(
              myrecv_indexes[index_of_jpid_in_myrecv_list],
              myrecv_counts[index_of_jpid_in_myrecv_list] * sizeof(int));
          myrecv_indexes[index_of_jpid_in_myrecv_list]
                        [myrecv_counts[index_of_jpid_in_myrecv_list] - 1] = j;
        }
        // do nothing if j is already included
      }
    }
  }

  // informing processors in my recv list about my needs
  for (int i = 0; i < myrecv_list_size; i++) {
    // send recv count
    MPI_Send(&myrecv_counts[i], 1, MPI_INT, myrecv_list[i], 0, MPI_COMM_WORLD);
    // send recv indexes
    MPI_Send(myrecv_indexes[i], myrecv_counts[i], MPI_INT, myrecv_list[i], 1,
             MPI_COMM_WORLD);
  }

  // informing processors not in my recv list that I need nothing from them
  for (int i = 0; i < numprocs; i++) {
    if (find(myrecv_list, myrecv_list_size, i) == -1) {
      int count = 0;
      MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }

  // getting informed about other processor's needs and constructing send list
  int mysend_list_size = 0;
  int *mysend_list = NULL;
  int *mysend_counts = NULL;
  int **mysend_indexes = NULL;
  for (int i = 0; i < numprocs; i++) {
    int count;
    MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (count != 0) { // processor i needs data from me
      mysend_list_size++;

      mysend_list = realloc(mysend_list, mysend_list_size * sizeof(int));
      mysend_list[mysend_list_size - 1] = i;

      mysend_counts = realloc(mysend_counts, myrecv_list_size * sizeof(int));
      mysend_counts[mysend_list_size - 1] = count;

      mysend_indexes =
          realloc(mysend_indexes, myrecv_list_size * sizeof(int *));
      mysend_indexes[mysend_list_size - 1] = malloc(count * sizeof(int));

      MPI_Recv(mysend_indexes[mysend_list_size - 1], count, MPI_INT, i, 1,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // constructing receive and send buffers
  double **my_recv_buffer = malloc(myrecv_list_size * sizeof(double *));
  for (int i = 0; i < myrecv_list_size; i++) {
    // initialize buffers to 0
    my_recv_buffer[i] = calloc(myrecv_counts[i], sizeof(double));
  }

  double **my_send_buffer = malloc(mysend_list_size * sizeof(double *));
  for (int i = 0; i < mysend_list_size; i++) {
    // initialize buffers to 0
    my_send_buffer[i] = calloc(mysend_counts[i], sizeof(double));
  }

  // create local X's
  double *local_X_new = malloc((LOCAL_ROW_COUNT) * sizeof(double));
  double *local_X_old = calloc((LOCAL_ROW_COUNT), sizeof(double));

  // execute the algorithm
  double error = DBL_MAX;
  while (error > EPSILON) {
    /* perform parallel SpMV */

    // exchange local status of X
    MPI_Request *reqs = malloc(myrecv_list_size * sizeof(MPI_Request));
    for (int i = 0; i < myrecv_list_size; i++) {
      MPI_Irecv(my_recv_buffer[i], myrecv_counts[i], MPI_DOUBLE, myrecv_list[i],
                0, MPI_COMM_WORLD, &reqs[i]);
    }

    for (int i = 0; i < mysend_list_size; i++) {
      MPI_Send(my_send_buffer[i], mysend_counts[i], MPI_DOUBLE, mysend_list[i],
               0, MPI_COMM_WORLD);
    }

    MPI_Waitall(myrecv_list_size, reqs, MPI_STATUS_IGNORE);

    free(reqs);

    for (int i = 0; i < LOCAL_ROW_COUNT; i++) {
      int displ = local_rows[0];
      double sum = 0;
      for (int k = local_rows[i]; k < local_rows[i + 1]; k++) {
        int j = local_columns[k - displ];
        // check if j is my
        int jpid = j / (LOCAL_ROW_COUNT); // mapping
        if (jpid != myid) {               // it is in my recvbuffer
          int index_of_jpid = find(myrecv_list, myrecv_list_size, jpid);
          int index_of_j = find(myrecv_indexes[index_of_jpid],
                                myrecv_counts[index_of_jpid], j);
          sum += local_values[k - displ] *
                 my_recv_buffer[index_of_jpid][index_of_j];
        } else {
          sum += local_values[k - displ] *
                 local_X_old[j - (myid * (LOCAL_ROW_COUNT))];
        }
      }
      local_X_new[i] = sum;
    }

    // remaining calculations
    double local_error_sum = 0;
    for (int i = 0; i < LOCAL_ROW_COUNT; i++) {
      local_X_new[i] = local_X_new[i] + local_b[i];
      local_error_sum += pow(local_X_new[i] - local_X_old[i], 2);
    }

    double total_error_sum;
    MPI_Allreduce(&local_error_sum, &total_error_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    error = sqrt(total_error_sum);

    // adjust my sendbuffer with updated values
    for (int i = 0; i < mysend_list_size; i++) {
      for (int j = 0; j < mysend_counts[i]; j++) {
        my_send_buffer[i][j] =
            local_X_new[mysend_indexes[i][j] - (myid * (LOCAL_ROW_COUNT))];
      }
    }

    // swap pointers
    double *temp = local_X_old;
    local_X_old = local_X_new;
    local_X_new = temp;
  }

  double *result_X;
  if (myid == 0) {
    result_X = malloc(matrix_size * sizeof(double));
  }

  MPI_Gather(local_X_old, LOCAL_ROW_COUNT, MPI_DOUBLE, result_X,
             LOCAL_ROW_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  double end = MPI_Wtime();
  double time = end - start;
  printf("Time: %f seconds\n", time);

  if (myid == 0) {
    char *output_vector_X_filename = argv[3];
    printResults(output_vector_X_filename, result_X, matrix_size);
  }

  free(local_rows);
  free(local_columns);
  free(local_values);
  free(local_b);

  free(myrecv_list);
  free(myrecv_counts);
  for (int i = 0; i < myrecv_list_size; i++) {
    free(myrecv_indexes[i]);
    free(my_recv_buffer[i]);
  }
  free(myrecv_indexes);
  free(my_recv_buffer);

  free(mysend_list);
  free(mysend_counts);
  for (int i = 0; i < mysend_list_size; i++) {
    free(mysend_indexes[i]);
    free(my_send_buffer[i]);
  }
  free(mysend_indexes);
  free(my_send_buffer);

  free(local_X_old);
  free(local_X_new);

  if (myid == 0) {
    free(result_X);
  }

  MPI_Finalize();

  return 0;
}
