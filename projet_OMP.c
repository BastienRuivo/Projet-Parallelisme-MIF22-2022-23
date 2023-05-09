#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <mpi.h>
#include <unistd.h>

void* mat_alloc(int M, int N)
{
  return calloc(M*N, sizeof(double));
}

/* A(i,j)<- 0 */
void mat_zero(int M, int N, double A[M][N])
{
  for (int i =0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      A[i][j] = 0;
    }
  }
}

/* A(i,j)<- 1 upper the diagonal */
void mat_shift(int M, int N, double A[M][N])
{
  for (int i =0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      A[i][j] = j == (i+1)%N ? 1 : 0;
    }
  }
}

/* A(i,j)<- random */
void mat_rand(int M, int N, double A[M][N])
{
  for (int i =0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      A[i][j] = drand48();
    }
  }
}

/* 1 iff A==B, else 0 */
int mat_is_eq(int M, int N, double A[M][N], double B[M][N])
{
  for (int i =0; i < M; ++i) {
    for (int j = 0; j < N; ++j) {
      if (A[i][j] != B[i][j]) return 0;
    }
  }
  return 1;
}


void mat_print(int M, int N, double A[M][N])
{
  for (int i =0; i < M; i++) {
    for (int j = 0; j < N; j++) {
      printf("%f ", A[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}


void kernel(int nLine, double in_mat[nLine][nLine], double out_mat[nLine][nLine], double conv_mask[3][3])
{
  for (int mat_x_idx = 1; mat_x_idx < nLine - 1; ++mat_x_idx) {
    for (int mat_y_idx = 1; mat_y_idx < nLine - 1; ++mat_y_idx) {
          out_mat[mat_x_idx][mat_y_idx] +=  
           ((in_mat[mat_x_idx - 1][mat_y_idx - 1] * conv_mask[0][0]) +
            (in_mat[mat_x_idx - 1][mat_y_idx] * conv_mask[0][1]) +
            (in_mat[mat_x_idx - 1][mat_y_idx + 1] * conv_mask[0][2]) +

            (in_mat[mat_x_idx][mat_y_idx - 1] * conv_mask[1][0]) +
            (in_mat[mat_x_idx][mat_y_idx] * conv_mask[1][1]) +
            (in_mat[mat_x_idx][mat_y_idx + 1] * conv_mask[1][2]) +

            (in_mat[mat_x_idx + 1][mat_y_idx - 1] * conv_mask[2][0]) +
            (in_mat[mat_x_idx + 1][mat_y_idx] * conv_mask[2][1]) +
            (in_mat[mat_x_idx + 1][mat_y_idx + 1] * conv_mask[2][2])) * 25724.0;
    }
  }
}



void kernel_MPI(int nLine, int nCol, double in_mat[nLine][nCol], double out_mat[nLine][nCol], double conv_mask[3][3])
{
  #pragma omp parallel for schedule(static) shared(in_mat, out_mat, conv_mask)
  for (int mat_x_idx = 1; mat_x_idx < nLine - 1; ++mat_x_idx) {
    for (int mat_y_idx = 1; mat_y_idx < nCol - 1; ++mat_y_idx) {
      out_mat[mat_x_idx][mat_y_idx] +=  
        ((in_mat[mat_x_idx - 1][mat_y_idx - 1] * conv_mask[0][0]) +
        (in_mat[mat_x_idx - 1][mat_y_idx] * conv_mask[0][1]) +  
        (in_mat[mat_x_idx - 1][mat_y_idx + 1] * conv_mask[0][2]) +

        (in_mat[mat_x_idx][mat_y_idx - 1] * conv_mask[1][0]) +
        (in_mat[mat_x_idx][mat_y_idx] * conv_mask[1][1]) +
        (in_mat[mat_x_idx][mat_y_idx + 1] * conv_mask[1][2]) +

        (in_mat[mat_x_idx + 1][mat_y_idx - 1] * conv_mask[2][0]) +
        (in_mat[mat_x_idx + 1][mat_y_idx] * conv_mask[2][1]) +
        (in_mat[mat_x_idx + 1][mat_y_idx + 1] * conv_mask[2][2])) * 25724.0;
    }
  }
}
/* Global variable for communication COMM_WORLD */
int my_rank = -1;
int num_process = -1;


int main(int argc, char const *argv[])
{
  /* MPI initialization */
  char hostname[MPI_MAX_PROCESSOR_NAME];
  int hostname_len;

  

  MPI_Init(0, 0);

  MPI_Get_processor_name(hostname, &hostname_len);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &num_process);

  if (argc < 2) {
    printf("usage : %s <SizeDomaine>\n", argv[0]);
    exit(1);
  }
  int N = atoi(argv[1]);
  int nLine = (N - 2) / num_process;
  if(nLine == 0) {
    printf("Size of the domain is too small for the number of process\n");
    exit(1);
  }
  int restLine = (N - 2) % num_process;
  int offset = 0;
  int addLine = 0;
  if(my_rank < restLine) {
    offset = my_rank;
    addLine = 1;
  } else {
    offset = restLine;
  }

  double (*in_mat)[N] = in_mat = mat_alloc(nLine + 2 + addLine, N);;
  double (*out_mat)[N] = mat_alloc(nLine + 2 + addLine, N);
  double (*in_complete)[N] = mat_alloc(N, N);
    struct timeval tval_before, tval_after, tval_result;

  if(my_rank == 0) {
    printf("Size of the domain : %d\n", N);
    printf("Number of process : %d\n", num_process);
    printf("Line size (it's N - 2 / num_process) %d\n", nLine);

    mat_rand(N, N, in_complete);
    gettimeofday(&tval_before, NULL);
    //#pragma omp parallel for
    for (size_t id = 1; id < num_process; id++)
    {
      int coffset = id < restLine ? id : restLine;
      if(id < restLine)
        MPI_Send(in_complete[id * (nLine) + coffset], N * (nLine + 3), MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
      else
        MPI_Send(in_complete[id * (nLine) + coffset], N * (nLine + 2), MPI_DOUBLE, id, 0, MPI_COMM_WORLD);
    }
    
    in_mat = in_complete;
    
  } else {
    MPI_Recv(in_mat, N * (nLine + 2 + addLine), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  double conv_mask [3][3] =
  {
    {-1, 0, 1},
    {-2, 0, 2},
    {-1, 0, 1}
  };

  kernel_MPI(nLine + 2 + addLine, N, in_mat, out_mat, conv_mask);

  if(my_rank != 0) {
    MPI_Send(out_mat[1], N * (nLine + addLine), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  else {
    double (* outm)[N] = mat_alloc(N, N);
    //#pragma omp parallel for
    for (size_t id = 1; id < num_process; id++)
    {
      int coffset = id < restLine ? id : restLine;
      if (id < restLine)
        MPI_Recv(outm[id * nLine + coffset + 1], N * (nLine + 1), MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      else
        MPI_Recv(outm[id * nLine + coffset + 1], N * nLine, MPI_DOUBLE, id, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    }

    # pragma omp parallel for collapse(2)
    for (size_t i = 0; i < nLine + addLine; i++)
    {
      for (size_t j = 0; j < N; j++)
      {
        outm[i + 1][j] = out_mat[i + 1][j];
      }
    }

    gettimeofday(&tval_after, NULL);
    timersub(&tval_after, &tval_before, &tval_result);
    
    double (* outm2)[N] = mat_alloc(N, N);
    kernel(N, in_complete, outm2, conv_mask);

    if (mat_is_eq(N, N, outm, outm2)) {
      printf("Success\n");
      printf("Time elapsed: %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);
    } else {
      printf("Failure\n");
    }
  }
  MPI_Finalize();
  return 0;
}

