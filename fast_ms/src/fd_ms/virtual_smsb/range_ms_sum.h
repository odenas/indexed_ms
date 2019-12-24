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
			uint64_t last_bit_idx,uint64_t * bitvec);
void range_ms_sum_fast_ext(uint32_t in_virtual_ms, uint32_t first_bit_idx,
			uint32_t last_bit_idx,uint32_t * bitvec,
			uint32_t * out_virtual_ms, uint32_t * ms_range_sum);

void range_ms_sum_fast_ext64(uint64_t in_virtual_ms, uint64_t first_bit_idx,
			uint64_t last_bit_idx,uint64_t * bitvec,
			uint64_t * out_virtual_ms, uint64_t * ms_range_sum);

static inline void range_ms_sum_fast_ext_word(uint32_t in_virtual_ms,
			uint32_t bitvec,uint32_t * out_virtual_ms,
			uint32_t * ms_range_sum)
{
	uint8_t i;
	uint8_t curr_byte;
	uint32_t range_sum=0;
	uint32_t virtual_ms=in_virtual_ms;
	range_ms_params_t range_ms_params;
	for(i=0;i<4;i++)
	{
		curr_byte=bitvec&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
		bitvec>>=8;
	};
	(*ms_range_sum)=range_sum;
	(*out_virtual_ms)=virtual_ms;
};

static inline void range_ms_sum_fast_ext_word64(uint64_t in_virtual_ms,
			uint64_t bitvec,uint64_t * out_virtual_ms,
			uint64_t * ms_range_sum)
{
	uint8_t i;
	uint8_t curr_byte;
	uint64_t range_sum=0;
	uint64_t virtual_ms=in_virtual_ms;
	range_ms_params_t range_ms_params;
	for(i=0;i<8;i++)
	{
		curr_byte=bitvec&0xff;
		range_ms_params=range_ms_lookup_table[curr_byte];
		range_sum+=virtual_ms*range_ms_params.sum_mul+range_ms_params.sum_add;
		virtual_ms+=range_ms_params.ms_incr;
		bitvec>>=8;
	};
	(*ms_range_sum)=range_sum;
	(*out_virtual_ms)=virtual_ms;
};


#endif
