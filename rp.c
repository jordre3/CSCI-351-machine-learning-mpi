/*
Copyright (c) 2016-2020 Jeremy Iverson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

/* assert */
#include <assert.h>

/* fabs */
#include <math.h>

/* MPI API */
#include <mpi.h>

#include <omp.h>

/* printf, fopen, fclose, fscanf, scanf */
#include <stdio.h>

/* EXIT_SUCCESS, malloc, calloc, free, qsort */
#include <stdlib.h>

#define MPI_SIZE_T MPI_UNSIGNED_LONG

struct distance_metric {
  size_t viewer_id;
  double distance;
};

static int
cmp(void const *ap, void const *bp)
{
  struct distance_metric const a = *(struct distance_metric*)ap;
  struct distance_metric const b = *(struct distance_metric*)bp;

  return a.distance < b.distance ? -1 : 1;
}

int
main(int argc, char * argv[])
{
  int ret, p, rank;
  size_t n, m, k;
  double * rating;

  /* Initialize MPI environment. */
  ret = MPI_Init(&argc, &argv);
  assert(MPI_SUCCESS == ret);

  /* Get size of world communicator. */
  ret = MPI_Comm_size(MPI_COMM_WORLD, &p);
  assert(ret == MPI_SUCCESS);

  /* Get my rank. */
  ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  assert(ret == MPI_SUCCESS);

  /* Validate command line arguments. */
  assert(2 == argc);

  /* Read input --- only if your rank 0. */
  if (0 == rank) {
    /* ... */
    char const * const fn = argv[1];

    /* Validate input. */
    assert(fn);

    /* Open file. */
    FILE * const fp = fopen(fn, "r");
    assert(fp);

    /* Read file. */
    fscanf(fp, "%zu %zu", &n, &m);

    rating = malloc(n * m * sizeof(*rating));

    /* Check for success. */
    assert(rating);

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < m; j++) {
        fscanf(fp, "%lf", &rating[i * m + j]);
      }
    }

    /* Close file. */
    ret = fclose(fp);
    assert(!ret);
  }
file:///usr/share/doc/HTML/en-US/index.html

  MPI_Bcast(&n, 1, MPI_SIZE_T, 0, MPI_COMM_WORLD);
  MPI_Bcast(&m, 1, MPI_SIZE_T, 0, MPI_COMM_WORLD);

  /* Compute base number of viewers. */
  size_t const base = 1 + ((n - 1) / p); // ceil(n / p)

  /* Compute local number of viewers. */
  size_t const ln = (rank + 1) * base > n ? n - rank * base : base;

  int * size;
  int * offset;

  /* Send viewer data to rest of processes. */
  if (0 == rank) {
    size = malloc(p * sizeof(int));
    assert(size);
    offset = malloc(p * sizeof(int));
    assert(offset);
    for (int r = 0; r < p; r++) {
      int const rn = (r + 1) * base > n ? n - r * base : base;
      size[r] = rn * m;
      offset[r] = r * base * m;
      // ret = MPI_Send(rating + r * base * m, rn * m, MPI_DOUBLE, r, 0,
      //   MPI_COMM_WORLD);
      // assert(MPI_SUCCESS == ret);
    }
  } else {
      rating = malloc(ln * m * sizeof(*rating));
      /* Check for success. */
      assert(rating);
    }
    ret = MPI_Scatterv(rating, size, offset, MPI_DOUBLE, rating, ln * m, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    assert(ret == MPI_SUCCESS);
  if (0 == rank){
    free(size);
    free(offset);
  }
  /* Allocate more memory. */
  double * const urating = malloc((m - 1) * sizeof(*urating));

  /* Check for success. */
  assert(urating);

  /* Get user input and send it to rest of processes. */
  if (0 == rank) {
    for (size_t j = 0; j < m - 1; j++) {
      printf("Enter your rating for movie %zu: ", j + 1);
      fflush(stdout);
      scanf("%lf", &urating[j]);
    }
  }

  MPI_Bcast(urating, m-1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  /* Allocate more memory. */
  //
  double * const distance = calloc(ln, sizeof(*distance));

  /* Check for success. */
  assert(distance);

  double ts = omp_get_wtime();
  /* Compute distances. */
  for (size_t i = 0; i < ln; i++) {
    for (size_t j = 0; j < m - 1; j++) {
      distance[i] += fabs(urating[j] - rating[i * m + j]);
    }
  }

  double te = omp_get_wtime();
  double local_elapsed = te - ts;
  double global_elapsed;
  ret = MPI_Reduce(&local_elapsed, &global_elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  assert(ret == MPI_SUCCESS);

  if(rank == 0){
    printf("Time: %lf\n", global_elapsed);
  }

  double * gatherDistances;
  int * recvCounts;
  int * displs;
  if(rank == 0){
    gatherDistances = malloc(n * sizeof(double));
    recvCounts = malloc(p * sizeof(int));
    displs = malloc(p * sizeof(int));
    for(size_t i = 0; i < p; i++){
      int const recvLen = (i + 1) * base > n ? n - i * base : base;
      recvCounts[i] = recvLen;
      displs[i] = i * base;
    }
   }

   ret = MPI_Gatherv(distance, ln, MPI_DOUBLE, gatherDistances, recvCounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   assert(ret == MPI_SUCCESS);

if(0 == rank){
  struct distance_metric * const finalDistances = calloc(n, sizeof(*finalDistances));
  for(int i = 0; i < n; i++){
    finalDistances[i].distance = gatherDistances[i];
    finalDistances[i].viewer_id = i;
  }
 qsort(finalDistances, n, sizeof(*finalDistances), cmp);

 /* Get user input. */
 printf("Enter the number of similar viewers to report: \n");
 scanf("%zu", &k);

 /* Output k viewers who are least different from the user. */
 printf("Viewer ID   Movie five   Distance\n");
 printf("---------------------------------\n");

 for (size_t i = 0; i < k; i++) {
   printf("%9zu   %10.1lf   %8.1lf\n", finalDistances[i].viewer_id + 1,
     rating[finalDistances[i].viewer_id * m + 4], finalDistances[i].distance);
 }

 printf("---------------------------------\n");

 /* Compute the average to make the prediction. */
 double sum = 0.0;
 for (size_t i = 0; i < k; i++) {
   sum += rating[finalDistances[i].viewer_id * m + 4];
 }

 /* Output prediction. */
 printf("The predicted rating for movie five is %.1lf.\n", sum / k);
 free(finalDistances);
}

  /* Deallocate memory. */
  free(rating);
  free(urating);
  free(distance);

  ret = MPI_Finalize();
  assert(MPI_SUCCESS == ret);

  return EXIT_SUCCESS;
}
