#ifndef naive_ms_range_max_h
#define naive_ms_range_max_h
#include<stdint.h>
#include"naive_ms.h"

uint32_t naive_ms_range_max(uint32_t first,uint32_t last,uint32_t bitvec_size,uint32_t * bitvec, uint32_t * _pos);
uint64_t naive_ms_range_max64(uint64_t first,uint64_t last,uint64_t bitvec_size,uint64_t * bitvec, uint64_t *_pos);
int8_t naive_ms_range_max_byte_incr(uint8_t byte,int8_t *_pos);

#endif
