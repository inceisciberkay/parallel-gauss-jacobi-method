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
  printf("Usage: ./gaussJacobi-mpi-p2 <input_matrix_A_csc> <input_vector_B> "
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

  int matrix_size, num_values, *columns, *rows;
  double *values;
  double *b;

  if (myid == 0) { // master process
    // reading input files
    char *input_matrix_A_filename = argv[1];
    char *input_vector_B_filename = argv[2];

    // read input matrix A
    if (!read_input_matrix_A_file(input_matrix_A_filename, &matrix_size,
                                  &num_values, &columns, &rows, &values)) {
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

  const int LOCAL_COLUMN_COUNT = matrix_size / numprocs;

  /* distribute input matrix A and input vector b to child processes */

  // partitioning columns, rows and values (uniform block column-wise
  // partitioning)
  int *send_counts_columns, *displs_columns, *send_counts_rows, *displs_rows,
      *send_counts_b, *displs_b;
  if (myid == 0) {
    send_counts_columns = malloc(numprocs * sizeof(int));
    displs_columns = malloc(numprocs * sizeof(int));

    send_counts_rows = malloc(numprocs * sizeof(int));
    displs_rows = malloc(numprocs * sizeof(int));

    send_counts_b = malloc(numprocs * sizeof(int));
    displs_b = malloc(numprocs * sizeof(int));

    for (int i = 0; i < numprocs; i++) {
      send_counts_columns[i] = LOCAL_COLUMN_COUNT + 1;
      displs_columns[i] = i * LOCAL_COLUMN_COUNT;

      send_counts_rows[i] = columns[(i + 1) * (LOCAL_COLUMN_COUNT)] -
                            columns[i * (LOCAL_COLUMN_COUNT)];
      displs_rows[i] =
          i == 0 ? 0 : displs_rows[i - 1] + send_counts_rows[i - 1];

      send_counts_b[i] = LOCAL_COLUMN_COUNT;
      displs_b[i] = i * LOCAL_COLUMN_COUNT;
    }
  }

  int rcvsize_columns = LOCAL_COLUMN_COUNT + 1;
  int *local_columns = malloc(rcvsize_columns * sizeof(int));
  MPI_Scatterv(columns, send_counts_columns, displs_columns, MPI_INT,
               local_columns, rcvsize_columns, MPI_INT, 0, MPI_COMM_WORLD);

  int rcvsize_rows = local_columns[rcvsize_columns - 1] - local_columns[0];
  int *local_rows = malloc(rcvsize_rows * sizeof(int));
  double *local_values = malloc(rcvsize_rows * sizeof(double));
  MPI_Scatterv(rows, send_counts_rows, displs_rows, MPI_INT, local_rows,
               rcvsize_rows, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatterv(values, send_counts_rows, displs_rows, MPI_DOUBLE, local_values,
               rcvsize_rows, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  // partitioning input vector b
  int rcvsize_b = LOCAL_COLUMN_COUNT;
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

  // constructing send list
  int my_send_list_size = 0;
  int *my_send_list = NULL;
  int *my_send_counts = NULL;
  int **my_send_indexes = NULL;

  for (int j = 0; j < LOCAL_COLUMN_COUNT; j++) {
    int displ = local_columns[0];
    for (int k = local_columns[j]; k < local_columns[j + 1]; k++) {
      int i = local_rows[k - displ];
      int ipid = i / (LOCAL_COLUMN_COUNT); // mapping
      if (ipid != myid) {
        int index_of_ipid_in_my_send_list =
            find(my_send_list, my_send_list_size, ipid);

        if (index_of_ipid_in_my_send_list == -1) {
          my_send_list_size++;

          my_send_list = realloc(my_send_list, my_send_list_size * sizeof(int));
          my_send_counts =
              realloc(my_send_counts, my_send_list_size * sizeof(int));
          my_send_indexes =
              realloc(my_send_indexes, my_send_list_size * sizeof(int *));

          index_of_ipid_in_my_send_list = my_send_list_size - 1;
          my_send_list[index_of_ipid_in_my_send_list] = ipid;
          my_send_counts[index_of_ipid_in_my_send_list] = 0;
          my_send_indexes[index_of_ipid_in_my_send_list] = NULL;
        }

        // check if i is already included in myrecv_indexes
        if (find(my_send_indexes[index_of_ipid_in_my_send_list],
                 my_send_counts[index_of_ipid_in_my_send_list], i) == -1) {
          // if not included
          my_send_counts[index_of_ipid_in_my_send_list]++;
          my_send_indexes[index_of_ipid_in_my_send_list] = realloc(
              my_send_indexes[index_of_ipid_in_my_send_list],
              my_send_counts[index_of_ipid_in_my_send_list] * sizeof(int));
          my_send_indexes[index_of_ipid_in_my_send_list]
                         [my_send_counts[index_of_ipid_in_my_send_list] - 1] =
                             i;
        }
        // do nothing if i is already included
      }
    }
  }

  // informing processors in my send list that I am going to send them data;
  // they must allocate receive buffers accordingly
  for (int i = 0; i < my_send_list_size; i++) {
    // send send count
    MPI_Send(&my_send_counts[i], 1, MPI_INT, my_send_list[i], 0,
             MPI_COMM_WORLD);
    // send send indexes
    MPI_Send(my_send_indexes[i], my_send_counts[i], MPI_INT, my_send_list[i], 1,
             MPI_COMM_WORLD);
  }

  // informing processors not in my send list that I will send nothing to them
  for (int i = 0; i < numprocs; i++) {
    if (find(my_send_list, my_send_list_size, i) == -1) {
      int count = 0;
      MPI_Send(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
    }
  }

  // getting informed about what other processors will send me and constructing
  // receive list accordingly
  int my_recv_list_size = 0;
  int *my_recv_list = NULL;
  int *my_recv_counts = NULL;
  int **my_recv_indexes = NULL;
  for (int i = 0; i < numprocs; i++) {
    int count;
    MPI_Recv(&count, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if (count != 0) { // processor i will send data to me
      my_recv_list_size++;

      my_recv_list = realloc(my_recv_list, my_recv_list_size * sizeof(int));
      my_recv_list[my_recv_list_size - 1] = i;

      my_recv_counts = realloc(my_recv_counts, my_recv_list_size * sizeof(int));
      my_recv_counts[my_recv_list_size - 1] = count;

      my_recv_indexes =
          realloc(my_recv_indexes, my_recv_list_size * sizeof(int *));
      my_recv_indexes[my_recv_list_size - 1] = malloc(count * sizeof(int));

      MPI_Recv(my_recv_indexes[my_recv_list_size - 1], count, MPI_INT, i, 1,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
  }

  // constructing receive and send buffers
  double **my_recv_buffer = malloc(my_recv_list_size * sizeof(double *));
  for (int i = 0; i < my_recv_list_size; i++) {
    // initialize buffers to 0
    my_recv_buffer[i] = calloc(my_recv_counts[i], sizeof(double));
  }

  double **my_send_buffer = malloc(my_send_list_size * sizeof(double *));
  for (int i = 0; i < my_send_list_size; i++) {
    // initialize buffers to 0
    my_send_buffer[i] = calloc(my_send_counts[i], sizeof(double));
  }

  // create local X's
  double *local_X_new = malloc((LOCAL_COLUMN_COUNT) * sizeof(double));
  double *local_X_old = calloc((LOCAL_COLUMN_COUNT), sizeof(double));

  // execute the algorithm
  double error = DBL_MAX;
  while (error > EPSILON) {
    /* perform parallel SpMV */

    // clear local X values
    for (int i = 0; i < LOCAL_COLUMN_COUNT; i++) {
      local_X_new[i] = 0;
    }

    // clear send buffer
    for (int i = 0; i < my_send_list_size; i++) {
      for (int j = 0; j < my_send_counts[i]; j++) {
        my_send_buffer[i][j] = 0.0;
      }
    }

    for (int j = 0; j < LOCAL_COLUMN_COUNT; j++) {
      int displ = local_columns[0];
      for (int k = local_columns[j]; k < local_columns[j + 1]; k++) {
        int i = local_rows[k - displ];
        // check if i is my; if so i need to write into my send buffer
        int ipid = i / LOCAL_COLUMN_COUNT; // mapping
        if (ipid != myid) {                // it is in my recvbuffer
          int index_of_ipid = find(my_send_list, my_send_list_size, ipid);
          int index_of_i = find(my_send_indexes[index_of_ipid],
                                my_send_counts[index_of_ipid], i);
          my_send_buffer[index_of_ipid][index_of_i] +=
              local_values[k - displ] * local_X_old[j];
        } else {
          local_X_new[i - myid * LOCAL_COLUMN_COUNT] +=
              local_values[k - displ] * local_X_old[j];
        }
      }
    }

    // exchange local status of X
    MPI_Request *reqs = malloc(my_recv_list_size * sizeof(MPI_Request));
    for (int i = 0; i < my_recv_list_size; i++) {
      MPI_Irecv(my_recv_buffer[i], my_recv_counts[i], MPI_DOUBLE,
                my_recv_list[i], 0, MPI_COMM_WORLD, &reqs[i]);
    }

    for (int i = 0; i < my_send_list_size; i++) {
      MPI_Send(my_send_buffer[i], my_send_counts[i], MPI_DOUBLE,
               my_send_list[i], 0, MPI_COMM_WORLD);
    }

    MPI_Waitall(my_recv_list_size, reqs, MPI_STATUS_IGNORE);

    free(reqs);

    // update my local X according to the values coming from other processes
    for (int i = 0; i < my_recv_list_size; i++) {
      for (int j = 0; j < my_recv_counts[i]; j++) {
        local_X_new[my_recv_indexes[i][j] - (myid) * (LOCAL_COLUMN_COUNT)] +=
            my_recv_buffer[i][j];
      }
    }

    // remaining calculations
    double local_error_sum = 0;
    for (int i = 0; i < LOCAL_COLUMN_COUNT; i++) {
      local_X_new[i] = local_X_new[i] + local_b[i];
      local_error_sum += pow(local_X_new[i] - local_X_old[i], 2);
    }

    double total_error_sum;
    MPI_Allreduce(&local_error_sum, &total_error_sum, 1, MPI_DOUBLE, MPI_SUM,
                  MPI_COMM_WORLD);

    error = sqrt(total_error_sum);

    // swap pointers
    double *temp = local_X_old;
    local_X_old = local_X_new;
    local_X_new = temp;
  }

  double *result_X;
  if (myid == 0) {
    result_X = malloc(matrix_size * sizeof(double));
  }

  MPI_Gather(local_X_old, LOCAL_COLUMN_COUNT, MPI_DOUBLE, result_X,
             LOCAL_COLUMN_COUNT, MPI_DOUBLE, 0, MPI_COMM_WORLD);

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

  free(my_recv_list);
  free(my_recv_counts);
  for (int i = 0; i < my_recv_list_size; i++) {
    free(my_recv_indexes[i]);
    free(my_recv_buffer[i]);
  }
  free(my_recv_indexes);
  free(my_recv_buffer);

  free(my_send_list);
  free(my_send_counts);
  for (int i = 0; i < my_send_list_size; i++) {
    free(my_send_indexes[i]);
    free(my_send_buffer[i]);
  }
  free(my_send_indexes);
  free(my_send_buffer);

  free(local_X_old);
  free(local_X_new);

  if (myid == 0) {
    free(result_X);
  }

  MPI_Finalize();

  return 0;
}
