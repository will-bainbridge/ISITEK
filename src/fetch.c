////////////////////////////////////////////////////////////////////////////////

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "fetch.h"

#define FETCH_MAX_STRING_LENGTH 1024

////////////////////////////////////////////////////////////////////////////////

struct s_FETCH
{
	int n_lines, max_n_lines; //number of data lines
	char *format; //data formats
	void **data; //the data
};

////////////////////////////////////////////////////////////////////////////////

FETCH fetch_new(char *format, int max_n_lines)
{
        //counters
        int i, j;

	//allocate the structure
	FETCH fetch = (FETCH)malloc(sizeof(struct s_FETCH));
	if(fetch == NULL) return NULL;

	//set the numbers
	fetch->n_lines = 0;
	fetch->max_n_lines = max_n_lines;

	//allocate and copy over the format string
	fetch->format = (char *)malloc((strlen(format) + 1) * sizeof(char));
	if(fetch->format == NULL) return NULL;
	strcpy(fetch->format,format);

	//size is the sum of all the data types
	int size = 0;
	for(i = 0; i < strlen(fetch->format); i ++)
	{
		switch(fetch->format[i])
		{
			case 'i': size += sizeof(int); break;
			case 'f': size += sizeof(float); break;
			case 'd': size += sizeof(double); break;
			case 'c': size += sizeof(char); break;
			case 's': size += sizeof(char*); break;
		}
	}

	//allocate the data
	fetch->data = (void **)malloc(fetch->max_n_lines * sizeof(void *));
	if(fetch->data == NULL) return NULL;
	fetch->data[0] = (void *)malloc(fetch->max_n_lines * size);
	if(fetch->data == NULL) return NULL;

	//pointer to the data
	void *d = fetch->data[0];

	//step through the data
	for(i = 0; i < fetch->max_n_lines; i ++)
	{
		fetch->data[i] = d;

		for(j = 0; j < strlen(fetch->format); j ++)
		{
			switch(fetch->format[j])
			{
				case 'i': d = (int*)d + 1; break;
				case 'f': d = (float*)d + 1; break;
				case 'd': d = (double*)d + 1; break;
				case 'c': d = (char*)d + 1; break;
				case 's': *((char**)d) = (char*)malloc(FETCH_MAX_STRING_LENGTH * sizeof(char));
					  if(*((char**)d) == NULL) return NULL;
					  d = (char**)d + 1; break;
			}
		}
	}

	//return array
	return fetch;
}

////////////////////////////////////////////////////////////////////////////////

int fetch_read(FILE *file, char *label, FETCH fetch)
{
	//check the file
	if(file == NULL) { return FETCH_FILE_ERROR; }
	rewind(file);

	//counters
	int i, offset, n_pieces = strlen(fetch->format);

	//pointer to the current value
	void *d;

	//reset the number of lines
	fetch->n_lines = 0;

	//allocate temporary storage
	char *line, *line_label, *line_data;
	line = (char *)malloc(FETCH_MAX_STRING_LENGTH * sizeof(char));
	if(line == NULL) return FETCH_MEMORY_ERROR;
	line_label = (char *)malloc(FETCH_MAX_STRING_LENGTH * sizeof(char));
	if(line_label == NULL) return FETCH_MEMORY_ERROR;
	line_data = (char *)malloc(FETCH_MAX_STRING_LENGTH * sizeof(char));
	if(line_data == NULL) return FETCH_MEMORY_ERROR;

	//read each line in turn
	while(fgets(line, FETCH_MAX_STRING_LENGTH, file) != NULL)
	{
		// get the first string on the line
		if(sscanf(line, "%s", line_label) == 1)
		{
			//check against the specified label
			if(strcmp(line_label, label) == 0)
			{
				//offset to the start of the data
				offset = strlen(label) + 1;

				//point to the first value
				d = fetch->data[fetch->n_lines];

				//loop over the desired bits of data
				for(i = 0; i < n_pieces; i ++)
				{
					//eat up whitespace
					while(line[offset] == ' ') offset ++;

					//read the data as a string
					if(sscanf(&line[offset], "%s", line_data) == 1)
					{
						//convert the string data to the desired type and increment the value pointer
						if(fetch->format[i] == 'i') {
							if(sscanf(line_data, "%i", (int*)d) != 1) break;
							d = (int*)d + 1;
						} else if(fetch->format[i] == 'f') {
							if(sscanf(line_data, "%f", (float*)d) != 1) break;
							d = (float*)d + 1;
						} else if(fetch->format[i] == 'd') {
							if(sscanf(line_data, "%lf", (double*)d) != 1) break;
							d = (double*)d + 1;
						} else if(fetch->format[i] == 'c') {
							if(sscanf(line_data, "%c", (char*)d) != 1) break;
							d = (char*)d + 1;
						} else if(fetch->format[i] == 's') {
							if(sscanf(line_data, "%s", *((char**)d)) != 1) break;
							d = (char**)d + 1;
						}

						//offset to the start of the next piece of data
						offset += strlen(line_data) + 1;
					}
				}

				//increment the number of lines if all values succesfully read
				if(i == n_pieces) (fetch->n_lines) ++;

				//quit if the maximum number of lines has been reached
				if(fetch->n_lines == fetch->max_n_lines) break;
			}
		}
	}

	//clean up and return the number of lines read
	free(line);
	free(line_label);
	free(line_data);
	return fetch->n_lines;
}

////////////////////////////////////////////////////////////////////////////////

void fetch_get(FETCH fetch, int line_index, int value_index, void *value)
{
	int i;
	void *d = fetch->data[line_index];

	for(i = 0; i < value_index; i ++)
	{
		switch(fetch->format[i])
		{
			case 'i': d = (int*)d + 1; break;
			case 'f': d = (float*)d + 1; break;
			case 'd': d = (double*)d + 1; break;
			case 'c': d = (char*)d + 1; break;
			case 's': d = (char**)d + 1; break;
		}
	}

	switch(fetch->format[value_index])
	{
		case 'i': *((int*)value) = *((int*)d); return;
		case 'f': *((float*)value) = *((float*)d); return;
		case 'd': *((double*)value) = *((double*)d); return;
		case 'c': *((char*)value) = *((char*)d); return;
		case 's': strcpy((char*)value,*((char**)d)); return;
	}
}

////////////////////////////////////////////////////////////////////////////////

void fetch_print(FETCH fetch)
{
	int i, j, n_pieces = strlen(fetch->format);
	void *d;

	for(i = 0; i < fetch->n_lines; i ++)
	{
		d = fetch->data[i];
		for(j = 0; j < n_pieces; j ++)
		{
			switch(fetch->format[j])
			{
				case 'i': printf("%i ",*((int*)d)); d = (int*)d + 1; break;
				case 'f': printf("%f ",*((float*)d)); d = (float*)d + 1; break;
				case 'd': printf("%lf ",*((double*)d)); d = (double*)d + 1; break;
				case 'c': printf("%c ",*((char*)d)); d = (char*)d + 1; break;
				case 's': printf("%s ",*((char**)d)); d = (char**)d + 1; break;
			}
		}
		printf("\n");
	}
}

////////////////////////////////////////////////////////////////////////////////

void fetch_destroy(FETCH fetch)
{
	int i, j, n_pieces = strlen(fetch->format);
	void *d;

	for(i = 0; i < fetch->max_n_lines; i ++)
	{
		d = fetch->data[i];
		for(j = 0; j < n_pieces; j ++)
		{
			switch(fetch->format[j])
			{
				case 'i': d = (int*)d + 1; break;
				case 'f': d = (float*)d + 1; break;
				case 'd': d = (double*)d + 1; break;
				case 'c': d = (char*)d + 1; break;
				case 's': free(*((char**)d)); 
					  d = (char**)d + 1; break;
			}
		}
	}

	free(fetch->format);
	free(fetch->data[0]);
	free(fetch->data);
	free(fetch);
}

////////////////////////////////////////////////////////////////////////////////

int fetch_value(FILE *file, char *label, char type, void *value)
{
	char format[2];

	format[0] = type;
	format[1] = '\0';

	FETCH fetch = fetch_new(format,1);
	if(fetch == NULL) return FETCH_MEMORY_ERROR;

	if(fetch_read(file, label, fetch) != 1) { fetch_destroy(fetch); return FETCH_FAIL; }

	fetch_get(fetch, 0, 0, value);

	fetch_destroy(fetch);

	return FETCH_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////

int fetch_vector(FILE *file, char *label, char type, int n, void *value)
{
	char *format = (char *)malloc((n + 1) * sizeof(char));
	if(format == NULL) return FETCH_MEMORY_ERROR;

	memset(format,type,n);
	format[n] = '\0';

	FETCH fetch = fetch_new(format,1);
	if(fetch == NULL) return FETCH_MEMORY_ERROR;

	if(fetch_read(file,label,fetch) != 1) { free(format); fetch_destroy(fetch); return FETCH_FAIL; }

	int i;
	void *d = fetch->data[0];
	void *v = value;

	for(i = 0; i < n; i ++)
	{
		switch(type)
		{
			case 'i': *((int*)v) = *((int*)d); v = (int*)v + 1; d = (int*)d + 1; break;
			case 'f': *((float*)v) = *((float*)d); v = (float*)v + 1; d = (float*)d + 1; break;
			case 'd': *((double*)v) = *((double*)d); v = (double*)v + 1; d = (double*)d + 1; break;
			case 'c': *((char*)v) = *((char*)d); v = (char*)v + 1; d = (char*)d + 1; break;
			case 's': strcpy(*((char**)v),*((char**)d)); v = (char**)v + 1; d = (char**)d + 1; break;
		}
	}

	free(format);
	fetch_destroy(fetch);

	return FETCH_SUCCESS;
}

////////////////////////////////////////////////////////////////////////////////
