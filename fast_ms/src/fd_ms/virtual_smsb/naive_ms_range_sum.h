#ifndef naive_ms_range_sum_h
#define naive_ms_range_sum_h
#include<stdint.h>
#include "naive_ms.h"


uint32_t naive_range_ms(uint32_t first,uint32_t last,uint32_t bitvec_size,uint32_t * bitvec);
uint64_t naive_range_ms64(uint64_t first,uint64_t last,uint64_t bitvec_size,const uint64_t * bitvec);

int32_t naive_range_ms_byte(int32_t virt_last_ms,uint8_t byte);

int8_t naive_ms_incr_byte(uint8_t byte);

#endif
