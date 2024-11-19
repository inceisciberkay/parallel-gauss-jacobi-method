#ifndef UTIL_H_
#define UTIL_H_

/*
 *		reads a a graph file in CSC format, each integer number
 *separated by line, First (numRows+1) lines define offsets array, the next
 *(numColumns) lines define edges array, and the next (numColumns) lines define
 *the weights array.
 *
 *
 *       This is used for testing on regular pc
 *
 *       @params:
 *		file_name: name of the file
 *		numRows: number of rows(vertices) in the file
 *		numColumns: number of columns (edges) in the file
 *      offsets: offsets(vertices) array
 *
 *       Returns 0 if there is a file format or a reading problem, 1 if
 *successful
 *
 *       Usage Example:
 *       int numVertices,numEdges,*offsets,*edges, *weights;
 *       int success =
 *read_file(argv[1],&numVertices,&numEdges,&offsets,&edges); if(success){ edges
 *= &(offsets[numVertices + 1]); weights = &(offsets[numVertices + 1 +
 *numEdges]);
 *           //Do the computation
 *       }
 **/
int read_input_matrix_A_file(char *file_name, int *numRows, int *numColumns,
                             int **vertices, int **edges, double **weights);

int read_input_vector_b_file(char *file_name, double **b);

void printResults(char *fileName, double *output_vector_X, int numVertices);

#endif /* UTIL_H_ */
