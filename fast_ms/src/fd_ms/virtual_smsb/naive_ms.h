#ifndef naive_ms_h
#define naive_ms_h
#include<stdint.h>
static inline uint32_t read_bit(uint32_t bit_idx,const uint32_t * bitvec)
{
	return (bitvec[bit_idx/32]>>(bit_idx%32))&1;
};

static inline uint64_t read_bit64(uint64_t bit_idx,const uint64_t * bitvec)
{
	return (bitvec[bit_idx/64]>>(bit_idx%64))&1;
};

uint32_t naive_select(uint32_t bit_idx,uint32_t bitvec_size,const uint32_t * bitvec);
uint64_t naive_select64(uint64_t bit_idx,uint64_t bitvec_size,const uint64_t * bitvec);

static inline uint32_t ms_at_idx(uint32_t bit_idx,uint32_t bitvec_size,uint32_t * bitvec)
{
	return naive_select(bit_idx,bitvec_size,bitvec)-2*bit_idx;
};

static inline uint64_t ms_at_idx64(uint64_t bit_idx,uint64_t bitvec_size,const uint64_t * bitvec)
{
	return naive_select64(bit_idx,bitvec_size,bitvec)-2*bit_idx;
};

int8_t naive_count1_byte(uint8_t byte);
int8_t naive_ms_incr_byte(uint8_t byte);

#endif
