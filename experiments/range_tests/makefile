
IID       = mm
DUMP_RIDX = ../../fast_ms/bin/dump_range_index.x
RANGE_Q   = ../../fast_ms/bin/range_queries_profile.x -ms_path $(IID).t.ms -from_max_idx 80000000
#../../fast_ms/bin/range_queries_profile.x -ms_path m.t.ms -ridx_path m.t.ms.4.range -block_size 0 -niter 20 -from_max_idx 100000000 -range_size 50000000



RANGE_SIZES = $(addsuffix 0000000,1 2 3 4 5 6 7 8 9 10)
BLOCK_SIZES = 0 4 8 16 32 64 128 256 512 1024 2048 4096
BLOCK_SIZES = 0 4 128 512 1024 4096

RANGES   = $(addsuffix .range,$(addprefix $(IID).t.ms.,$(BLOCK_SIZES)))
TARGETS  = $(subst .range,.csv,$(RANGES))

.SECONDARY : $(RANGES)

all :  $(TARGETS)

$(TARGETS) : $(IID).t.ms.%.csv : $(IID).t.ms.%.range
	echo "block_size,range_size,nqueries,method,time_ms" >$@
	for i in $(RANGE_SIZES); do \
		$(RANGE_Q) -ridx_path $< -range_size $$i -block_size $(subst $(IID).t.ms.,,$(subst .csv,,$@)) -niter 10 >> $@;\
	done

$(IID).t.ms.%.range : $(IID).t.ms
	$(DUMP_RIDX) -ms_path $<  -block_size $(subst $(IID).t.ms.,,$(subst .range,,$@)) -time_usage 1 >>$(IID).t.ms.range.csv

clean : 
	rm -f  $(TARGETS) $(IID).t.ms.range.csv
