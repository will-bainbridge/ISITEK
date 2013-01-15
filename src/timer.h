//////////////////////////////////////////////////////////////////

#include <time.h>

//////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////
