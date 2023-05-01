#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>


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

void kernel(int N, double in_mat[N][N], double out_mat[N][N], double conv_mask[3][3]);


int main(int argc, char const *argv[])
{

  if (argc < 2) {
    printf("usage : %s <SizeDomaine>\n", argv[0]);
    exit(1);
  }
  int N = atoi(argv[1]);

  double (*in_mat)[N] = mat_alloc(N, N);
  double (*out_mat)[N] = mat_alloc(N, N);

  mat_rand(N, N, in_mat);
  mat_zero(N, N, out_mat);

  double conv_mask [3][3] =
  {
    {-1, 0, 1},
    {-2, 0, 2},
    {-1, 0, 1}
  };

  // Call kernel, measure time
  //
  struct timeval tval_before, tval_after, tval_result;
  gettimeofday(&tval_before, NULL);

  kernel(N, in_mat, out_mat, conv_mask);

  gettimeofday(&tval_after, NULL);
  timersub(&tval_after, &tval_before, &tval_result);

  //mat_print(N, N, in_mat);
  //mat_print(N, N, out_mat);
  printf("Time needed: %ld.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);


  // Check kernel correctness
  //
  double (*test_in_mat)[5] = mat_alloc(5, 5);
  double (*test_out_mat)[5] = mat_alloc(5, 5);
  double test_ok_mat [5][5] = 
  {
    {0, 0, 0, 0, 0},
    {0, 51448, 0, -51448, 0},
    {0, 25724, 51448, 0, 0},
    {0, -25724, 25724,51448, 0},
    {0, 0, 0, 0, 0},
  };
  mat_zero(5, 5, test_in_mat);
  mat_shift(5, 5, test_in_mat);
  kernel(5, test_in_mat, test_out_mat, conv_mask);
  if (!mat_is_eq(5, 5, test_out_mat, test_ok_mat)) {
    fprintf(stderr,"*** incorrect result! ***\n");
    exit(1);
  }

  return 0;
}

void kernel(int N, double in_mat[N][N], double out_mat[N][N], double conv_mask[3][3])
{
  for (int mat_x_idx = 1; mat_x_idx < N-1; ++mat_x_idx) {
    for (int mat_y_idx = 1; mat_y_idx < N-1; ++mat_y_idx) {
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
