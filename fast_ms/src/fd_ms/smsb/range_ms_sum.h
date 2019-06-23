#ifndef range_ms_sum_h
#define range_ms_sum_h
#include<stdint.h>
typedef struct 
{
	uint8_t sum_mul;
	int8_t sum_add;
	int8_t ms_incr;
} range_ms_params_t;
extern range_ms_params_t range_ms_lookup_table[256];
uint32_t range_ms_sum_fast(uint32_t first_ms, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec);
uint64_t range_ms_sum_fast64(uint64_t first_ms, uint64_t first_bit_idx,
			uint64_t last_bit_idx, const uint64_t * bitvec);


#endif
