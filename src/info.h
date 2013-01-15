//////////////////////////////////////////////////////////////////

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

//////////////////////////////////////////////////////////////////
