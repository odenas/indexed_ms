#include<stdlib.h>
#include<stdio.h>
#include"naive_ms_range_sum.h"

typedef struct 
{
	uint8_t sum_mul;
	int8_t sum_add;
	int8_t ms_incr;
}range_ms_params_t;

range_ms_params_t range_ms_lookup_table[256];

int main()
{
	uint32_t byte_value;
	int32_t ms1,ms2;
	FILE *f0=fopen("ms_range_sum_tables.c","w");
	fprintf(f0,"typedef struct {uint8_t sum_mul;uint8_t sum_add;int8_t ms_incr;} range_ms_params_t;\n");
	fprintf(f0,"range_ms_params_t range_ms_lookup_table[256]={");
	for(byte_value=0;byte_value<=255;byte_value++)
	{
		ms1=naive_range_ms_byte(1,byte_value);
		ms2=naive_range_ms_byte(2,byte_value);
		range_ms_lookup_table[byte_value].sum_mul=ms2-ms1;
		range_ms_lookup_table[byte_value].sum_add=2*ms1-ms2;
		range_ms_lookup_table[byte_value].ms_incr=
			naive_ms_incr_byte(byte_value);
		if((byte_value%8)==0)
			fprintf(f0,"\n");
		fprintf(f0,"{%d,%d,%d}",range_ms_lookup_table[byte_value].sum_mul,
				range_ms_lookup_table[byte_value].sum_add,
				range_ms_lookup_table[byte_value].ms_incr);
		if(byte_value<255)
			fprintf(f0,",");
	};
	fprintf(f0,"};\n");
	fclose(f0);
	return 0;
};
