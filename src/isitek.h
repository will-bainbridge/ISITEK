#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expression.h"
#include "fetch.h"
#include "memory.h"
#include "sparse.h"
#include "structure.h"

// error handling
#define exit_if_false(value,message) \
	do { \
		if(!(value)) \
		{ \
			fprintf(stderr,"\n[ERROR %s:%i] %s\n\n",__FILE__,__LINE__,message); \
			exit(EXIT_FAILURE); \
		} \
	} while(0)

// warning handling
#define warn_if_false(value,message) \
	do { \
		if(!(value)) \
		{ \
			fprintf(stderr,"\n[WARNING %s:%i] %s\n",__FILE__,__LINE__,message); \
		} \
	} while(0)

// timing utility
#define print_time(format,function) \
	do { \
		struct timespec handle_time[2]; \
		clock_gettime(CLOCK_MONOTONIC, &handle_time[0]); \
		function; \
		clock_gettime(CLOCK_MONOTONIC, &handle_time[1]); \
		double handle_seconds = difftime(handle_time[1].tv_sec, handle_time[0].tv_sec); \
		long handle_nano_seconds = handle_time[1].tv_nsec - handle_time[0].tv_nsec; \
		fprintf(stdout,format, handle_seconds + handle_nano_seconds / 1.0e9); \
		fflush(stdout); \
	} while(0)

// maximum numbers for allocation purposes
#define MAX_FACE_N_NODES 2
#define MAX_FACE_N_BOUNDARIES 2
#define MAX_FACE_N_BORDERS 2
#define MAX_ELEMENT_N_FACES 4
#define MAX_STRING_LENGTH 1024
#define MAX_EXPRESSION_N_RECURSIONS 10

// expression ordering
#define EXPRESSION_LOCATION_INDEX 0
#define EXPRESSION_LOCATION_LABELS "xy"
#define EXPRESSION_NORMAL_INDEX 2
#define EXPRESSION_NORMAL_LABELS "ab"
#define EXPRESSION_VARIABLE_INDEX 4

// calculate number of basis functions from order
#define ORDER_TO_N_BASIS(x) ((x)*(x+1)/2)
#define ORDER_TO_N_GAUSS(x) (2*(x)-1)
#define ORDER_TO_N_HAMMER(x) 12 // (7+2*(x>2)+3*(x>3))

// utilities
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
