#include "mpi.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

//#define NUM_POINTS 524288 
#define NUM_POINTS 2048

unsigned int X_axis[NUM_POINTS];
unsigned int Y_axis[NUM_POINTS];

typedef struct {
	unsigned int * prim;
	unsigned int * secd;
} CoorArray;
 
int partition (unsigned int *, unsigned int *, int, int);
void local_quick_sort (unsigned int *, unsigned int *, int, int);
void print_local_chunk(unsigned int *, unsigned int *, int, int, int);
CoorArray local_merge(unsigned int *, unsigned int *,
		                  unsigned int *, unsigned int *, int);
CoorArray tree_merge(unsigned int *, unsigned int *, int, int, int);
 
void find_quadrants (num_quadrants)
     int num_quadrants;
{
  /* YOU NEED TO FILL IN HERE */








}

void parallel_sort (numprocs, myid)
	int numprocs;
	int myid;
{
	// calculate and broadcast the chunk size
	int size = (NUM_POINTS + numprocs - 1) / numprocs;
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// allocate memory and scatter the x and y coordiante
	unsigned int *chunk_X;
	unsigned int *chunk_Y;
	chunk_X = (unsigned int *)calloc(size, sizeof(unsigned int));
	chunk_Y = (unsigned int *)calloc(size, sizeof(unsigned int));
	MPI_Scatter(X_axis, size, MPI_INT, chunk_X, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(Y_axis, size, MPI_INT, chunk_Y, size, MPI_INT, 0, MPI_COMM_WORLD);

	
	local_quick_sort(chunk_X, chunk_Y, 0, size - 1);
	CoorArray res_sorted_x;
	res_sorted_x = tree_merge(chunk_X, chunk_Y, size, numprocs, myid);

	local_quick_sort(chunk_Y, chunk_X, 0, size - 1);
	CoorArray res_sorted_y;
	res_sorted_y = tree_merge(chunk_Y, chunk_X, size, numprocs, myid);

	/*
	if (myid == 0) {
		int i;
		fprintf(stdout, "Check result: \n");
		for (i = 0; i < NUM_POINTS - 1; i++) {
			if (res_sorted_x.prim[i] > res_sorted_x.prim[i + 1]) {
				fprintf(stdout, "Check failed @ %d: %d %d\n", i, res_sorted_x.prim[i], res_sorted_x.prim[i + 1]);
			}
			//fprintf(stdout, "%12d ", res.res_x[i]);
		}
		fprintf(stdout, "\n");
	}
	*/
	return;
}

void print_local_chunk(chunk_X, chunk_Y, size, numprocs, myid)
	unsigned int *chunk_X;
	unsigned int *chunk_Y;
	int size;
	int numprocs;
	int myid;
{
	int id;
	for (id = 0; id < numprocs; id++) {
		if (id == myid) {
			int i;
			fprintf (stdout, "process %d:\n", myid);
			for (i = 0; i < size; i++) {
				fprintf (stdout, "%d ", chunk_X[i]);
			}
			fprintf (stdout, "\n");
			for (i = 0; i < size; i++) {
				fprintf (stdout, "%d ", chunk_Y[i]);
			}
			fprintf (stdout, "\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
}
			

void local_quick_sort (prim_axis, secd_axis, left, right)
	unsigned int *prim_axis;
	unsigned int *secd_axis;
	int left;
	int right;
{
	int middle;

	if (left < right) {
		middle = partition(prim_axis, secd_axis, left, right);
		local_quick_sort(prim_axis, secd_axis, left, middle - 1);
		local_quick_sort(prim_axis, secd_axis, middle + 1, right);
	}

	return;
}

int partition (prim_axis, secd_axis, left, right)
	unsigned int *prim_axis;
	unsigned int *secd_axis;
	int left;
	int right;
{
	int i, j;
	unsigned int pivot, t;

	// since data are random numbers, fixed pivot is ok
	pivot = prim_axis[left];
	i = left;
	j = right + 1;

	while(1) {
		do ++i; while (prim_axis[i] <= pivot && i <= right);
		do --j; while (prim_axis[j] > pivot);
		if (i > j) break;
		t = prim_axis[i];
		prim_axis[i] = prim_axis[j];
		prim_axis[j] = t;
		t = secd_axis[i];
		secd_axis[i] = secd_axis[j];
		secd_axis[j] = t;
	}
	t = prim_axis[left];
	prim_axis[left] = prim_axis[j];
	prim_axis[j] = t;
	t = secd_axis[left];
	secd_axis[left] = secd_axis[j];
	secd_axis[j] = t;
	
	return j;
}

CoorArray local_merge(a1_x, a1_y, a2_x, a2_y, size)
	unsigned int * a1_x;
	unsigned int * a1_y;
	unsigned int * a2_x;
	unsigned int * a2_y;
	int size;
{
	unsigned int * res_x;
	unsigned int * res_y;
	res_x = (unsigned int *)malloc(size * 2 * sizeof(unsigned int));
	res_y = (unsigned int *)malloc(size * 2 * sizeof(unsigned int));

	int i = 0;
	int j = 0;
	int k = 0;
	while (i < size && j < size) {
		if (a1_x[i] < a2_x[j]) {
			res_x[k] = a1_x[i];
			res_y[k] = a1_y[i];
			i++;
			k++;
		} else {
			res_x[k] = a2_x[j];
			res_y[k] = a2_y[j];
			j++;
			k++;
		}
	}
	if (i == size) {
		while (j < size) {
			res_x[k] = a2_x[j];
			res_y[k] = a2_y[j];
			j++;
			k++;
		}
	} else {
		while (i < size) {
			res_x[k] = a1_x[i];
			res_y[k] = a1_y[i];
			i++;
			k++;
		}
	}

	CoorArray res;
	res.prim = res_x;
	res.secd = res_y;
	return res;
}

CoorArray tree_merge(a_x, a_y, size, numprocs, myid)
	 unsigned int * a_x;
	 unsigned int * a_y;
	 int size;
	 int numprocs;
	 int myid;
{
	int step = 1;
	unsigned int *t_x;
	unsigned int *t_y;

	while (step < numprocs) {
		int s = size * step;
		if (myid % (2 * step) == 0) {
			if (myid + step < numprocs) {
				t_x = (unsigned int *)malloc(s * sizeof(unsigned int));
				t_y = (unsigned int *)malloc(s * sizeof(unsigned int));
				MPI_Status status;
				MPI_Recv(t_x, s, MPI_INT, myid + step, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(t_y, s, MPI_INT, myid + step, 1, MPI_COMM_WORLD, &status);
				
				CoorArray result;
				result = local_merge(a_x, a_y, t_x, t_y, s);
				/*
				if (myid == 4) {
					int i;
					fprintf(stdout, "step: %d\n", step);
					for (i = 0; i < s * 2; i++) {
						fprintf(stdout, "%d ", result.res_x[i]);
					}
					fprintf(stdout, "\n");
				}
				*/
				if (step > 1) {
					free(a_x);
					free(a_y);
				}
				a_x = result.prim;
				a_y = result.secd;
				free(t_x);
				free(t_y);
			}
		} else {
				int target_id = myid - step;
				MPI_Send(a_x, s, MPI_INT, target_id, 0, MPI_COMM_WORLD);
				MPI_Send(a_y, s, MPI_INT, target_id, 1, MPI_COMM_WORLD);
				if (step > 1) {
					free(a_x);
					free(a_y);
				}
				break;
		}
		step *= 2;
		//MPI_Barrier(MPI_COMM_WORLD);
	}

	CoorArray res;
	res.prim = a_x;
	res.secd = a_y;
	return res;
}

int main(argc,argv)
  int argc;
 char *argv[];
{
  int num_quadrants;
  int myid, numprocs;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
    
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);

  if (argc != 2)
    {
      fprintf (stderr, "Usage: recursive_bisection <#of quadrants>\n");
      MPI_Finalize();
      exit (0);
    }

  fprintf (stderr,"Process %d on %s\n", myid, processor_name);

  num_quadrants = atoi (argv[1]);

  if (myid == 0)
    fprintf (stdout, "Extracting %d quadrants with %d processors \n", num_quadrants, numprocs);

  if (myid == 0)
    {
      int i;

      srand (10000);
      
      for (i = 0; i < NUM_POINTS; i++)
				X_axis[i] = (unsigned int)rand();

      for (i = 0; i < NUM_POINTS; i++)
				Y_axis[i] = (unsigned int)rand();
    }

  //MPI_Bcast(&X_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  
  //MPI_Bcast(&Y_axis, NUM_POINTS, MPI_INT, 0, MPI_COMM_WORLD);  

  //find_quadrants (num_quadrants);
	parallel_sort (numprocs, myid);
 
  MPI_Finalize();
  return 0;
}
  

