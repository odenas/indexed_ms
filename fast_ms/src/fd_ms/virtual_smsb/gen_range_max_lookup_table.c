#include<stdlib.h>
#include<stdio.h>
#include"naive_ms_range_max.h"

typedef struct 
{
	int8_t max_incr;
	int8_t max_pos;
	int8_t ms_incr;
	uint8_t count1;
}ms_range_max_params_t;

ms_range_max_params_t ms_range_max_lookup_table[256];

int main()
{
	uint32_t byte_value;
	int8_t max_pos;
	int8_t max_ms;
	FILE *f0=fopen("ms_range_max_tables.c","w");
	fprintf(f0,"typedef struct {int8_t max_incr;int8_t max_pos;int8_t ms_incr;uint8_t count1;} ms_range_max_params_t;\n");
	fprintf(f0,"ms_range_max_params_t ms_range_max_lookup_table[256]={");
	for(byte_value=0;byte_value<=255;byte_value++)
	{
		max_ms=naive_ms_range_max_byte_incr(byte_value,&max_pos);
		ms_range_max_lookup_table[byte_value].max_incr=max_ms;
		ms_range_max_lookup_table[byte_value].max_pos=max_pos;
		ms_range_max_lookup_table[byte_value].ms_incr=
			naive_ms_incr_byte(byte_value);
		ms_range_max_lookup_table[byte_value].count1=
			naive_count1_byte(byte_value);
		if((byte_value%8)==0)
			fprintf(f0,"\n");
		fprintf(f0,"{%d,%d,%d,%d}",ms_range_max_lookup_table[byte_value].max_incr,
				ms_range_max_lookup_table[byte_value].max_pos,
				ms_range_max_lookup_table[byte_value].ms_incr,
				ms_range_max_lookup_table[byte_value].count1);
		if(byte_value<255)
			fprintf(f0,",");
	};
	fprintf(f0,"};\n");
	fclose(f0);
	return 0;
};
