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

// printing
#define print_info(...) \
	do { \
		printf("\033[1;32m   INFO \033[0m"); \
		printf(__VA_ARGS__); \
		printf("\n"); \
	} while(0)
#define print_output(...) \
	do { \
		printf("\033[1;36m OUTPUT \033[0m"); \
		printf(__VA_ARGS__); \
		printf("\n"); \
	} while(0)
#define print_continue(...) \
	do { \
		printf("        "); \
		printf(__VA_ARGS__); \
		printf("\n"); \
	} while(0)

// error handling
#define exit_if_false(value,...) \
	do { \
		if(!(value)) \
		{ \
			printf("\033[1;31m  ERROR \033[0m\033[1m%s:%i \033[0m",__FILE__,__LINE__); \
			printf(__VA_ARGS__); \
			printf("\n"); \
			exit(EXIT_FAILURE); \
		} \
	} while(0)

// warning handling
#define warn_if_false(value,...) \
	do { \
		if(!(value)) \
		{ \
			printf("\033[1;33mWARNING \033[0m\033[1m%s:%i \033[0m",__FILE__,__LINE__); \
			printf(__VA_ARGS__); \
			printf("\n"); \
		} \
	} while(0)

// timing utilities
struct timespec timer_time[2];
#define timer_elapsed() (difftime(timer_time[1].tv_sec, timer_time[0].tv_sec) + (timer_time[1].tv_nsec - timer_time[0].tv_nsec) / 1.0e9)
#define timer_start() \
	do { \
		clock_gettime(CLOCK_MONOTONIC, &timer_time[0]); \
		printf("\033[1;35m   TIME \033[0m\033[1mSTART \033[0m\n"); \
	} while(0)
#define timer_reset() \
	do { \
		clock_gettime(CLOCK_MONOTONIC, &timer_time[1]); \
		printf("\033[1;35m   TIME \033[0m%.6e s \033[1mRESET \033[0m\n",timer_elapsed()); \
		timer_time[0] = timer_time[1]; \
	} while(0)
#define timer_print() \
	do { \
		clock_gettime(CLOCK_MONOTONIC, &timer_time[1]); \
		printf("\033[1;35m   TIME \033[0m%.6e s\n",timer_elapsed()); \
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
#define ORDER_TO_N_GAUSS(x) (x)
#define ORDER_TO_N_HAMMER(x) 12 // (7+2*(x>2)+3*(x>3))

// utilities
#define MIN(x,y) (x < y ? x : y)
#define MAX(x,y) (x > y ? x : y)
