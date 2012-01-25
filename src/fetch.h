#ifndef FETCH_H
#define FETCH_H

#define FETCH_FILE_ERROR -1
#define FETCH_MEMORY_ERROR -2

#define FETCH_SUCCESS 1
#define FETCH_FAIL 0

typedef struct s_FETCH * FETCH;

FETCH fetch_new(char *format, int max_n_lines);
int fetch_read(FILE *file, char *label, FETCH fetch);
void fetch_get(FETCH fetch, int line_index, int value_index, void *value);
void fetch_print(FETCH fetch);
void fetch_destroy(FETCH fetch);
int fetch_value(FILE *file, char *label, char type, void *value);
int fetch_vector(FILE *file, char *label, char type, int n, void *value);

#endif
