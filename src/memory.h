#ifndef MEMORY_H
#define MEMORY_H

#define MAX_STRING_LENGTH 1024

int ** matrix_integer_new(int **old, int height, int width);
double ** matrix_double_new(double **old, int height, int width);
char ** matrix_character_new(char **old, int height, int width);
void matrix_free(void **matrix);

int *** tensor_integer_new(int ***old, int height, int width, int depth);
double *** tensor_double_new(double ***old, int height, int width, int depth);
void tensor_free(void ***tensor);

#endif
