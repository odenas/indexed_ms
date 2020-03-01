#ifndef ms_range_max_h
#define ms_range_max_h
#include<stdint.h>
typedef struct
{
	int8_t max_incr;
	int8_t max_pos;
	int8_t ms_incr;
	uint8_t count1;
} ms_range_max_params_t;

extern ms_range_max_params_t ms_range_max_lookup_table[256];

uint32_t ms_range_max_fast(uint32_t first_ms_idx, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec, uint32_t * _pos);
uint64_t ms_range_max_fast64(uint64_t first_ms_idx, uint64_t first_bit_idx,
			uint64_t last_bit_idx,const uint64_t * bitvec, uint64_t * _pos);


#endif
